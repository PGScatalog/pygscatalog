"""Score parity with saved PLINK 2 and pgsc_calc v2.3.0 references.

v3 uses Python, DuckDB and Zarr for matching and scoring. A reference from v3
could repeat the same defect on both sides of the comparison. Standalone PLINK
checks allele-weight arithmetic independently. The v2 workflow also checks
compatibility with its matching, scoring and aggregation stages. Both references
use PLINK. They check different pipeline behaviours but do not provide two
independent scoring algorithms.

Committed VCFs and raw .sscore files make the comparison reproducible. The test
needs no network, containers or installed PLINK version. It combines the original
outputs instead of using a manually prepared total. PLINK and Nextflow run only
during deliberate reference regeneration. The test calls the current v3 stages
in the same process.

The synthetic cohort isolates allele counting. It uses complete calls, one
additive score, REF and ALT effects, signed weights, multiallelic SNPs and two
indels. It excludes missingness, imputation, PARs and heteroplasmy. These limits
help locate failures without combining different biological conventions.
See tests/data/v2_parity/README.md for the fixture and regeneration commands."""

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
    """Test LOAD through COMBINE with multiallelic calls and different ploidies."""
    is_v2 = reference_name == "pgsc_calc_v2.3.0"
    reference = DATA_DIR / reference_name
    # v2's pinned matcher excludes MT. Compare the same supported variants in both
    # versions. The standalone PLINK case still covers all chromosomes.
    chroms = CHROMS[:-1] if is_v2 else CHROMS
    n_matched = 50 if is_v2 else 52
    filename = "synthetic_{chrom}_additive_0.sscore" if is_v2 else "chr{chrom}.sscore"
    sample_column = "IID" if is_v2 else "#IID"
    sum_column = "PGS999997_SUM" if is_v2 else "WEIGHT_SUM"
    # Autosomes are diploid in every sample. Non-PAR X is diploid in XX samples and
    # haploid in XY samples. Y and MT calls are haploid. Both references use --xchr-
    # model 1. PLINK doubles male X dosage by default. v3 counts the single allele in
    # the VCF call.
    sample_ids = {"XX_1", "XX_2", "XY_1", "XY_2"}
    expected = defaultdict(lambda: defaultdict(Decimal))
    chromosome_expected = {}
    # Check the required files and samples before summing. A missing shard or
    # duplicate row could otherwise produce a plausible but incorrect total.
    assert {p.name for p in reference.glob("*.sscore")} == {
        filename.format(chrom=chrom) for chrom in chroms
    }
    for chrom in chroms:
        with (reference / filename.format(chrom=chrom)).open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        # XX samples have no Y contribution and need no Y imputation. The Y file
        # deliberately excludes them. The current reader cannot distinguish chromosome
        # absence from a missing haploid call with VCF '.'.
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
            # Dosage counts copies of the named effect allele, not all ALT alleles. At
            # a multiallelic 1/2 call, an ALT2 effect has dosage 1. A REF effect has
            # dosage 0. The weighted sum must use the correct allele. The unweighted
            # dosage sum helps detect errors that positive and negative weights could
            # cancel.
            for field, plink_field in (
                ("score", sum_column),
                ("dosage_sum", "NAMED_ALLELE_DOSAGE_SUM"),
            ):
                contribution = Decimal(row[plink_field])
                assert contribution.is_finite(), (chrom, key, field)
                chromosome_expected[chrom][key][field] = contribution
                expected[key][field] += contribution
    # The average denominator counts observed alleles, not variants or effect-allele
    # copies. All calls are present and neither run imputes. Thus, PLINK's ALLELE_CT
    # equals DENOM here. Missing calls or imputation could change that equality.
    # Chromosomes contribute different allele counts, so do not average WEIGHT_AVG
    # values. Divide the combined sums and round once.
    for values in expected.values():
        values["score_avg"] = values["score"] / values["allele_count"]
        for field in ("score", "dosage_sum", "score_avg"):
            values[field] = values[field].quantize(
                Decimal("0.000001"), rounding=ROUND_HALF_UP
            )

    if is_v2:
        # Agreement with PLINK shards alone could miss a v2 aggregation error. Also
        # check the workflow's original final output. This check ensures that
        # reconstructed totals represent the reported v2 results.
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

    # v3 reads ploidy from each GT call. The declared sex map checks those calls but
    # does not determine dosage. Two X variants cannot give a conclusive sex check.
    # This fixture tests scoring, not sex-check sensitivity.
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

    # run does not yet expose --keep_multiallelic. Use the public stages. Retain
    # multiallelic matches so the test cannot silently reduce coverage to biallelic
    # variants. MATCH processes each chromosome. MATCH_COMBINE applies the overlap
    # gate after it combines all chromosome matches.
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
        # The current gate uses >. Require all supported variants. For v2, 49/52 <
        # 0.95 < 50/52. For the full fixture, 51/52 < 0.99 < 52/52.
        "--min_overlap",
        "0.95" if is_v2 else "0.99",
        "--out_dir",
        combined,
    )
    if is_v2:
        # Equal totals can hide different allele choices. Compare selected variants
        # directly. v2's join uses only ID and duplicates multiallelic candidates. It
        # can report an incorrect ALT annotation. matched_effect_allele identifies the
        # allele sent to PLINK.
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

        # v3's log derives its allele string from the score row and match type. Check
        # the index that SCORE uses against the original VCF. A correct log must not
        # hide an incorrect ALT1/ALT2 selection.
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
        # Every sample has 98 alleles in the full fixture, or 96 without MT. XY
        # samples have two fewer X alleles and two additional Y alleles. Incorrect X
        # and Y counts could therefore cancel in the total. Score and dosage errors
        # can also cancel. Swapping chr2 and chr8 weights preserves totals for these
        # repeated genotypes. Check all three fields per chromosome, including the
        # absence of XX rows on Y.
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
    # A dictionary comparison can hide duplicate output rows. Check row count and
    # identity so each sample receives exactly one score.
    assert len(actual) == len(expected) == 4
    assert set(expected) == keys
    assert {(r["accession"], r["sample_id"]) for r in actual} == keys
    assert {r["sampleset"] for r in actual} == {SAMPLESET}
    assert all(int(r["n_matched"]) == n_matched for r in actual)
    # n_matched describes the whole score even though XX samples have no Y
    # contribution. allele_count and dosage_sum describe each sample. Weights are
    # multiples of 1/16, so binary arithmetic represents these sums exactly. Compare
    # the public six-decimal output with Decimal. A floating-point tolerance could
    # hide a counting error. Only division for score_avg requires rounding.
    for row in actual:
        key = (row["accession"], row["sample_id"])
        for field in ("score", "allele_count", "dosage_sum", "score_avg"):
            assert Decimal(row[field]) == expected[key][field], (key, field)
