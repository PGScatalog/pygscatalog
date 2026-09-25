"""Score parity with saved PLINK 2 and pgsc_calc v2.3.0 references.

The v3 rewrite moves matching and scoring into Python/DuckDB/Zarr. A reference
computed by v3 itself could preserve the same bug on both sides of a comparison.
Standalone PLINK checks allele-weight arithmetic independently; the actual v2
workflow additionally checks compatibility with the earlier matching, scoring
and aggregation pipeline. Both references use PLINK, so these are complementary
checks of pipeline behavior, not two independent scoring algorithms.

The committed VCFs and raw .sscore files make the comparison reproducible without
network access, container runtimes or changes to installed PLINK versions. This
test aggregates those original outputs rather than trusting a manually prepared
total. PLINK and Nextflow run only when deliberately regenerating the references;
the code under test is the current v3 pipeline, called in-process.

The small synthetic cohort isolates allele counting: fully called genotypes,
one additive score, REF and ALT effects, signed weights, multiallelic SNPs and
two indels. Missingness, imputation, PARs and heteroplasmy are intentionally absent
so failures can be traced without also resolving those biological conventions.
See tests/data/v2_parity/README.md for the fixture and regeneration commands.
"""

from __future__ import annotations

import csv
import gzip
import pathlib
from collections import defaultdict
from decimal import ROUND_HALF_UP, Decimal

import polars as pl
import pysam
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
    # v2's pinned matcher filters out MT. Compare the same supported variants on
    # both sides instead of treating its missing mitochondrial contribution as
    # a v3 regression; the standalone PLINK case still covers the complete set.
    chroms = CHROMS[:-1] if is_v2 else CHROMS
    n_matched = 50 if is_v2 else 52
    filename = "synthetic_{chrom}_additive_0.sscore" if is_v2 else "chr{chrom}.sscore"
    sample_column = "IID" if is_v2 else "#IID"
    sum_column = "PGS999997_SUM" if is_v2 else "WEIGHT_SUM"
    # In this fixture, autosomes are diploid in everyone, non-PAR X is diploid
    # in XX samples and haploid in XY samples, and Y/MT calls are haploid.
    # Both references explicitly use --xchr-model 1: PLINK's default doubles
    # male X dosage, whereas v3 counts the single allele recorded in the VCF.
    sample_ids = {"XX_1", "XX_2", "XY_1", "XY_2"}
    expected = defaultdict(lambda: defaultdict(Decimal))
    chromosome_expected = {}
    # Enumerate the required files and samples before summing: a missing shard
    # or duplicate row could otherwise become a plausible but incorrect total.
    assert {p.name for p in reference.glob("*.sscore")} == {
        filename.format(chrom=chrom) for chrom in chroms
    }
    for chrom in chroms:
        with (reference / filename.format(chrom=chrom)).open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        # XX samples have no Y contribution, not an unknown Y genotype to impute.
        # Their absence from the Y file is deliberate: VCF '.' cannot distinguish
        # chromosome absence from a missing haploid call in the current reader.
        wanted = {"XY_1", "XY_2"} if chrom == "Y" else sample_ids
        assert len(rows) == len(wanted), chrom
        assert {r[sample_column] for r in rows} == wanted, chrom
        chromosome_expected[chrom] = {}
        for row in rows:
            key = ("PGS999997", row[sample_column])
            count = int(row["ALLELE_CT"])
            assert count > 0, (chrom, key)
            chromosome_expected[chrom][key] = {"allele_count": count}
            expected[key]["allele_count"] += count
            # Dosage counts copies of the named effect allele, not all ALT
            # alleles. At a multiallelic 1/2 call, an ALT2 effect has dosage 1;
            # a REF effect has dosage 0. The weighted sum must respect that
            # identity, while the unweighted dosage sum helps expose mistakes
            # that positive and negative weights could otherwise cancel.
            for field, plink_field in (
                ("score", sum_column),
                ("dosage_sum", "NAMED_ALLELE_DOSAGE_SUM"),
            ):
                contribution = Decimal(row[plink_field])
                assert contribution.is_finite(), (chrom, key, field)
                chromosome_expected[chrom][key][field] = contribution
                expected[key][field] += contribution
    # The average's denominator is the number of observed alleles, not variants
    # or effect-allele copies. Because every call is present and neither run
    # imputes, PLINK's ALLELE_CT also equals its scoring DENOM here. Do not assume
    # that equivalence for a future fixture with missing calls or imputation.
    # Chromosomes contribute different allele counts, so averaging WEIGHT_AVG
    # would weight them incorrectly. Divide the combined sums and round once.
    for values in expected.values():
        values["score_avg"] = values["score"] / values["allele_count"]
        for field in ("score", "dosage_sum", "score_avg"):
            values[field] = values[field].quantize(
                Decimal("0.000001"), rounding=ROUND_HALF_UP
            )

    if is_v2:
        # Agreement with PLINK shards alone would miss a v2 aggregation issue.
        # Retain the workflow's original final output as a second check that
        # the reconstructed totals represent what v2 actually reported.
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

    # v3 observes ploidy from each GT call; the declared sex map validates those
    # calls rather than determining dosage. Two X variants cannot make its sex
    # check conclusive: this fixture tests scoring, not sex-check sensitivity.
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
    # Retaining multiallelic matches is essential: silently dropping them would
    # turn this into biallelic-only coverage. MATCH sees each chromosome; the
    # overlap gate is applied only after pooling all chromosome matches.
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
    if is_v2:
        # Equal totals can conceal different allele choices. Compare selected
        # variants directly, not v2's exploded candidate rows: its ID-only join
        # duplicates multiallelic candidates and can report the wrong ALT, while
        # matched_effect_allele is the allele actually passed to PLINK.
        with gzip.open(reference / "synthetic_log.csv.gz", "rt") as handle:
            selected_v2 = [
                r for r in csv.DictReader(handle) if r["match_status"] == "matched"
            ]
        selected_v3 = []
        for path in (combined / "logs").rglob("*.csv.gz"):
            with gzip.open(path, "rt") as handle:
                selected_v3.extend(
                    r for r in csv.DictReader(handle) if r["match_summary"] == "matched"
                )
        expected_matches = {
            (r["accession"], int(r["row_nr"])): (
                r["chr_name"],
                int(r["chr_position"]),
                r["matched_effect_allele"],
            )
            for r in selected_v2
        }
        actual_matches = {
            (r["accession"], int(r["row_nr"])): (
                r["chr_name"],
                int(r["chr_position"]),
                r["matched_effect_allele"],
            )
            for r in selected_v3
        }
        assert len(selected_v2) == len(expected_matches) == n_matched
        assert len(selected_v3) == len(actual_matches) == n_matched
        assert actual_matches == expected_matches

        # v3's log derives its allele string from the score row and match type.
        # Also resolve the index consumed by SCORE against the original VCF so
        # a correct-looking log cannot hide an incorrect ALT1/ALT2 selection.
        target_alleles = {}
        for path, chrom in targets:
            with pysam.VariantFile(path) as vcf:
                for variant in vcf:
                    target_alleles[chrom, variant.pos] = variant.alleles
        selected_table = [
            r
            for r in pl.read_parquet(combined / "match_table.parquet").to_dicts()
            if r["match_result"]["is_matched"]
        ]
        table_keys = {(r["accession"], r["row_nr"]) for r in selected_table}
        assert len(selected_table) == len(table_keys) == n_matched
        assert table_keys == set(expected_matches)
        for row in selected_table:
            key = (row["accession"], row["row_nr"])
            chrom, position, effect = expected_matches[key]
            index = row["match_result"]["effect_allele_idx"]
            alleles = target_alleles[chrom, position]
            assert 0 <= index < len(alleles), key
            assert alleles[index] == effect, key

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
        # Every sample happens to have 98 alleles in the full fixture (96 without
        # MT): XY samples' two fewer X alleles are offset by their two Y alleles.
        # A wrong X count and wrong Y count could therefore cancel in the total.
        # Scores and dosage sums can also conceal compensating errors: swapping
        # chr2 and chr8 weights preserves totals for these repeated genotypes.
        # Check all three fields per chromosome, including no XX rows on Y.
        partial_rows = pl.read_parquet(partials[-1]).to_dicts()
        assert len(partial_rows) == len(chromosome_expected[chrom]), chrom
        assert {(r["accession"], r["sample_id"]) for r in partial_rows} == set(
            chromosome_expected[chrom]
        ), chrom
        for row in partial_rows:
            key = (row["accession"], row["sample_id"])
            for field, value in chromosome_expected[chrom][key].items():
                assert Decimal(str(row[field])) == value, (chrom, key, field)
    stage("combine", "--partial_scores", *partials, "--out_dir", tmp_path / "combined")

    with gzip.open(tmp_path / "combined" / "scores.txt.gz", "rt") as handle:
        actual = list(csv.DictReader(handle))
    keys = {("PGS999997", sample) for sample in sample_ids}
    # A dictionary comparison alone can hide duplicate output rows. Check row
    # count and identity as well, so every sample receives exactly one score.
    assert len(actual) == len(expected) == 4
    assert set(expected) == keys
    assert {(r["accession"], r["sample_id"]) for r in actual} == keys
    assert {r["sampleset"] for r in actual} == {SAMPLESET}
    assert all(int(r["n_matched"]) == n_matched for r in actual)
    # n_matched is score-wide even though XX samples have no Y contribution;
    # allele_count and dosage_sum describe what each sample actually carries.
    # Weights are multiples of 1/16, so these sums are binary-exact. Compare the
    # public six-decimal output with Decimal instead of allowing a float tolerance
    # to hide a counting error; only the division for score_avg needs rounding.
    for row in actual:
        key = (row["accession"], row["sample_id"])
        for field in ("score", "allele_count", "dosage_sum", "score_avg"):
            assert Decimal(row[field]) == expected[key][field], (key, field)
