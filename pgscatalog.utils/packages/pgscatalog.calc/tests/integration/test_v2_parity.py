"""Score parity with saved PLINK 2 and pgsc_calc v2.3.0 references.

PLINK and Nextflow are only used to regenerate references, never in this test.
See tests/data/v2_parity/README.md for the fixture and regeneration commands.
"""

from __future__ import annotations

import csv
import gzip
import pathlib
from collections import defaultdict
from decimal import ROUND_HALF_UP, Decimal

import polars as pl
import pytest

from .conftest import SAMPLESET, load_targets, run_cli

DATA_DIR = pathlib.Path(__file__).parent.parent / "data" / "v2_parity"
CHROMS = (*map(str, range(1, 23)), "X", "Y", "MT")
SCORE_FILE = DATA_DIR / "score.txt.gz"


@pytest.mark.parametrize("reference_name", ["plink_reference", "pgsc_calc_v2.3.0"])
def test_chromosomes_match_reference_scores(reference_name, tmp_path, monkeypatch):
    """Exercise LOAD through COMBINE with multiallelic and mixed-ploidy calls."""
    is_v2 = reference_name == "pgsc_calc_v2.3.0"
    reference = DATA_DIR / reference_name
    # v2's pinned matcher supports 1-22, X and Y, but filters out MT.
    chroms = CHROMS[:-1] if is_v2 else CHROMS
    n_matched = 50 if is_v2 else 52
    filename = "synthetic_{chrom}_additive_0.sscore" if is_v2 else "chr{chrom}.sscore"
    sample_column = "IID" if is_v2 else "#IID"
    sum_column = "PGS999997_SUM" if is_v2 else "WEIGHT_SUM"
    sample_ids = {"XX_1", "XX_2", "XY_1", "XY_2"}
    expected = defaultdict(lambda: defaultdict(Decimal))
    chromosome_counts = {}
    assert {p.name for p in reference.glob("*.sscore")} == {
        filename.format(chrom=chrom) for chrom in chroms
    }
    for chrom in chroms:
        with (reference / filename.format(chrom=chrom)).open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        wanted = {"XY_1", "XY_2"} if chrom == "Y" else sample_ids
        assert len(rows) == len(wanted), chrom
        assert {r[sample_column] for r in rows} == wanted, chrom
        chromosome_counts[chrom] = {}
        for row in rows:
            key = ("PGS999997", row[sample_column])
            count = int(row["ALLELE_CT"])
            assert count > 0, (chrom, key)
            chromosome_counts[chrom][key] = count
            expected[key]["allele_count"] += count
            for field, plink_field in (
                ("score", sum_column),
                ("dosage_sum", "NAMED_ALLELE_DOSAGE_SUM"),
            ):
                contribution = Decimal(row[plink_field])
                assert contribution.is_finite(), (chrom, key, field)
                expected[key][field] += contribution
    # Round once after adding chromosome sums, never average WEIGHT_AVG.
    for values in expected.values():
        values["score_avg"] = values["score"] / values["allele_count"]
        for field in ("score", "dosage_sum", "score_avg"):
            values[field] = values[field].quantize(
                Decimal("0.000001"), rounding=ROUND_HALF_UP
            )

    if is_v2:
        # The saved end-to-end v2 output must also agree with its raw shards.
        with gzip.open(reference / "aggregated_scores.txt.gz", "rt") as handle:
            aggregate = list(csv.DictReader(handle, delimiter="\t"))
        assert len(aggregate) == 4
        assert {(r["PGS"], r["IID"]) for r in aggregate} == set(expected)
        for row in aggregate:
            key = (row["PGS"], row["IID"])
            assert row["sampleset"] == "synthetic"
            for field, column in (
                ("score", "SUM"),
                ("allele_count", "DENOM"),
                ("score_avg", "AVG"),
            ):
                value = Decimal(row[column]).quantize(
                    Decimal("0.000001"), rounding=ROUND_HALF_UP
                )
                assert value == expected[key][field], (key, column)

    targets = [(DATA_DIR / f"target.chr{chrom}.vcf.gz", chrom) for chrom in chroms]
    zips = load_targets(
        targets,
        tmp_path / "load",
        SCORE_FILE,
        monkeypatch,
        sex_map=DATA_DIR / "sex.tsv",
    )

    def stage(command, *args):
        run_cli(["pgsc_calc", command, *map(str, args)], monkeypatch)

    # run does not expose --keep_multiallelic yet; exercise the public stages.
    matches, samples, partials = [], [], []
    for chrom, archive in zip(chroms, zips, strict=True):
        out = tmp_path / "match" / chrom
        stage(
            "match",
            "--zarr_zip_file",
            archive,
            "--score_paths",
            SCORE_FILE,
            "--keep_multiallelic",
            "--out_dir",
            out,
        )
        matches.append(out / "matches.parquet")
        samples.append(out / "samples.parquet")
    combined = tmp_path / "matched"
    stage(
        "match-combine",
        "--matches",
        *matches,
        "--samples",
        *samples,
        "--score_paths",
        SCORE_FILE,
        "--sampleset",
        SAMPLESET,
        # Current gate uses >. Require all supported variants: 49/52 < 0.95
        # < 50/52 for v2, and 51/52 < 0.99 < 52/52 for the full fixture.
        "--min_overlap",
        "0.95" if is_v2 else "0.99",
        "--out_dir",
        combined,
    )
    for chrom, archive in zip(chroms, zips, strict=True):
        out = tmp_path / "score" / chrom
        stage(
            "score",
            "--zarr_zip_file",
            archive,
            "--score_paths",
            SCORE_FILE,
            "--match_table",
            combined / "match_table.parquet",
            "--out_dir",
            out,
        )
        partials.append(out / "partial_scores.parquet")
        # X and Y counts can cancel in the combined denominator: check each
        # chromosome independently, including the absence of XX samples on Y.
        partial_rows = pl.read_parquet(partials[-1]).to_dicts()
        assert len(partial_rows) == len(chromosome_counts[chrom]), chrom
        assert {
            (r["accession"], r["sample_id"]): r["allele_count"] for r in partial_rows
        } == chromosome_counts[chrom], chrom
    stage("combine", "--partial_scores", *partials, "--out_dir", tmp_path / "combined")

    with gzip.open(tmp_path / "combined" / "scores.txt.gz", "rt") as handle:
        actual = list(csv.DictReader(handle))
    keys = {("PGS999997", sample) for sample in sample_ids}
    assert len(actual) == len(expected) == 4
    assert set(expected) == keys
    assert {(r["accession"], r["sample_id"]) for r in actual} == keys
    assert {r["sampleset"] for r in actual} == {SAMPLESET}
    assert all(int(r["n_matched"]) == n_matched for r in actual)
    # Check every numerical output against PLINK, with no float tolerance.
    for row in actual:
        key = (row["accession"], row["sample_id"])
        for field in ("score", "allele_count", "dosage_sum", "score_avg"):
            assert Decimal(row[field]) == expected[key][field], (key, field)
