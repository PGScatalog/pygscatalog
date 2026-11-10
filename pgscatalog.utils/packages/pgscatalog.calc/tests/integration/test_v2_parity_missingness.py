"""Check the differences in v2 and v3 policies for missing calls.

With at least 50 samples, v2 and v3 fill missing calls from observed allele
frequencies (ADR 0009). They intentionally use different denominators for the
average. v3 counts observed alleles. A sample with no calls has an imputed sum
but no defined average. v2 also counts imputed alleles.

Donor calls give integer dosages of 0, 1 and 2. The fixture has no DS/GP inputs
or fractional dosage examples. The cohort meets the actual minimum of 50 samples."""

import csv
import gzip
import runpy
from decimal import Decimal
from pathlib import Path

import polars as pl
import pytest

from .parity_helpers import read_reference, run_case, validate_inputs

DATA = Path(__file__).parent.parent / "data/v2_parity/missingness"
ACCESSION = "PGS990010"


def test_missingness_is_a_documented_policy_difference(tmp_path, monkeypatch):
    validate_inputs(DATA)
    reference = read_reference(DATA / "pgsc_calc_v2.3.0")
    validate_inputs(DATA, DATA / "pgsc_calc_v2.3.0/no_mean_imputation")
    diagnostic = read_reference(DATA / "pgsc_calc_v2.3.0/no_mean_imputation")
    targets = [(DATA / f"target.chr{chrom}.vcf.gz", chrom) for chrom in ("1", "2")]
    actual = run_case(targets, [DATA / "score.txt.gz"], tmp_path, monkeypatch)
    samples = {f"S{i:02d}" for i in range(52)}
    assert len(actual) == 52
    assert {(r["accession"], r["sample_id"]) for r in actual.to_dicts()} == {
        (ACCESSION, sample) for sample in samples
    }
    assert set(reference) == {
        (chrom, ACCESSION, sample) for chrom in ("1", "2") for sample in samples
    }
    assert set(diagnostic) == set(reference)
    with gzip.open(DATA / "pgsc_calc_v2.3.0/aggregated_scores.txt.gz", "rt") as handle:
        aggregate = list(csv.DictReader(handle, delimiter="\t"))
    assert len(aggregate) == 52
    assert {(r["PGS"], r["IID"]) for r in aggregate} == {
        (ACCESSION, sample) for sample in samples
    }
    for row in aggregate:
        assert Decimal(row["SUM"]) == Decimal("2.5")
        assert Decimal(row["DENOM"]) == 8
        assert Decimal(row["AVG"]) == Decimal("0.3125")
    with gzip.open(
        DATA / "pgsc_calc_v2.3.0/no_mean_imputation/aggregated_scores.txt.gz", "rt"
    ) as handle:
        diagnostic_aggregate = list(csv.DictReader(handle, delimiter="\t"))
    assert len(diagnostic_aggregate) == 52
    assert {(r["PGS"], r["IID"]) for r in diagnostic_aggregate} == {
        (ACCESSION, sample) for sample in samples
    }
    for row in diagnostic_aggregate:
        shards = [diagnostic[chrom, ACCESSION, row["IID"]] for chrom in ("1", "2")]
        score = sum(r["score"] for r in shards)
        denom = sum(r["denom"] for r in shards)
        assert Decimal(row["SUM"]) == score
        assert Decimal(row["DENOM"]) == denom
        if denom:
            assert Decimal(row["AVG"]) == score / denom
        else:
            assert row["AVG"] == ""

    # Heterozygous donors give one copy at the first locus. Homozygous REF donors give
    # zero ALT copies on chr1 and two REF copies on chr2. These cases detect an
    # incorrect imputer that always fills missing calls with one copy.
    for index, (chrom, score, dosage, missing) in enumerate(
        (
            ("1", Decimal("2"), 1, {"S00", "S01"}),
            ("2", Decimal("0.5"), 3, {"S00", "S02"}),
        )
    ):
        partial = pl.read_parquet(
            tmp_path / "score" / str(index) / "partial_scores.parquet"
        )
        assert len(partial) == 52
        assert set(partial["sample_id"]) == samples
        for row in partial.to_dicts():
            sample = row["sample_id"]
            observed = int(sample not in missing)
            raw = reference[chrom, ACCESSION, sample]
            assert raw["score"] == score
            # PLINK reports observed named alleles even when its weighted sum includes
            # imputation. v3 includes imputed calls in its dosage statistic.
            assert raw["dosage_sum"] == dosage * observed
            assert raw["allele_count"] == 4 * observed
            assert raw["denom"] == 4
            no_mean = diagnostic[chrom, ACCESSION, sample]
            assert no_mean["score"] == score * observed
            assert no_mean["dosage_sum"] == dosage * observed
            assert no_mean["allele_count"] == no_mean["denom"] == 4 * observed
            assert row["accession"] == ACCESSION
            assert Decimal(str(row["score"])) == score
            assert row["dosage_sum"] == dosage
            assert row["allele_count"] == 4 * observed

    for row in actual.to_dicts():
        count = (
            0
            if row["sample_id"] == "S00"
            else 4
            if row["sample_id"] in {"S01", "S02"}
            else 8
        )
        assert row["n_matched"] == 4
        assert Decimal(str(row["score"])) == Decimal("2.5")
        assert row["allele_count"] == count
        assert row["dosage_sum"] == 4
        if count:
            assert Decimal(str(row["score_avg"])) == Decimal("2.5") / count
        else:
            assert row["score_avg"] is None

    selected = pl.read_parquet(tmp_path / "matched/match_table.parquet").to_dicts()
    assert len(selected) == 4
    assert {
        (r["row_nr"], r["match_result"]["effect_allele_idx"])
        for r in selected
        if r["match_result"]["is_matched"]
    } == {(0, 1), (1, 1), (2, 0), (3, 0)}
    with gzip.open(DATA / "pgsc_calc_v2.3.0/synthetic_log.csv.gz", "rt") as handle:
        matched = [r for r in csv.DictReader(handle) if r["match_status"] == "matched"]
    assert len(matched) == 4
    assert {(int(r["row_nr"]), r["matched_effect_allele"]) for r in matched} == {
        (0, "C"),
        (1, "C"),
        (2, "A"),
        (3, "A"),
    }


@pytest.mark.parametrize(
    "filename,error",
    [
        ("all_uncalled.chr1.vcf.gz", "no called alleles across the sampleset"),
        ("too_small.chr1.vcf.gz", "n_samples=4 is less than n_minimum_samples=50"),
    ],
)
def test_imputation_rejects_unsupported_inputs(filename, error, tmp_path, monkeypatch):
    # Entirely uncalled sites must survive LOAD and fail imputation explicitly.
    # Four samples with missing calls fail the minimum sample count.
    # Neither case may produce fabricated zero scores.
    with pytest.raises(ValueError, match=error):
        run_case(
            [(DATA / filename, "1")], [DATA / "score.txt.gz"], tmp_path, monkeypatch
        )
    assert not (tmp_path / "score/0/partial_scores.parquet").exists()
    assert not (tmp_path / "combined/scores.txt.gz").exists()
    if filename.startswith("all_uncalled"):
        assert not pl.read_parquet(tmp_path / "match/0/matches.parquet").is_empty()


@pytest.mark.parametrize("case,n_files", [("missingness", 11), ("zero_weight", 7)])
def test_missingness_and_zero_weight_inputs_regenerate(case, n_files, tmp_path):
    source = DATA.parent / case
    # Regeneration must never replace the independently computed references.
    reference = tmp_path / "pgsc_calc_v2.3.0"
    reference.mkdir()
    marker = reference / "external.sscore"
    marker.write_bytes(b"external reference\n")
    runpy.run_path(str(source / "generate.py"))["generate"](tmp_path)
    files = [path for path in tmp_path.iterdir() if path.is_file()]
    assert len(files) == n_files
    for path in files:
        assert path.read_bytes() == (source / path.name).read_bytes(), path.name
    assert marker.read_bytes() == b"external reference\n"
