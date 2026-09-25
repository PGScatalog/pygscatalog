"""A zero coefficient does not remove a matched locus from score statistics.

Allele count describes observed alleles at matched scoring rows, independently
of effect magnitude. The same applies to dosage sum. An all-zero score therefore
still has a defined zero average for called samples, rather than a null average.
"""

import csv
import gzip
from decimal import Decimal
from pathlib import Path

import polars as pl

from .parity_helpers import read_reference, run_case, validate_inputs

DATA = Path(__file__).parent.parent / "data/v2_parity/zero_weight"


def test_zero_weights_preserve_matched_allele_statistics(tmp_path, monkeypatch):
    validate_inputs(DATA)
    reference = read_reference(DATA / "pgsc_calc_v2.3.0")
    accessions = ("PGS990040", "PGS990041")
    actual = run_case(
        [(DATA / "target.chr1.vcf.gz", "1")],
        [DATA / f"{accession}.txt.gz" for accession in accessions],
        tmp_path,
        monkeypatch,
    )
    expected = {
        (accession, f"S{i}"): Decimal(2 * (2 - i) if accession == "PGS990040" else 0)
        for accession in accessions
        for i in range(3)
    }
    assert set(reference) == {("1", *key) for key in expected}
    partial = pl.read_parquet(tmp_path / "score/0/partial_scores.parquet")
    for frame in (partial, actual):
        assert len(frame) == 6
        assert {(r["accession"], r["sample_id"]) for r in frame.to_dicts()} == set(
            expected
        )
        for row in frame.to_dicts():
            key = row["accession"], row["sample_id"]
            assert Decimal(str(row["score"])) == expected[key]
            assert row["allele_count"] == 4, key
            assert row["dosage_sum"] == 2, key
            raw = reference["1", *key]
            assert raw["score"] == expected[key]
            assert raw["allele_count"] == raw["denom"] == 4
            assert raw["dosage_sum"] == 2
    for row in actual.to_dicts():
        key = row["accession"], row["sample_id"]
        assert row["n_matched"] == 2
        assert Decimal(str(row["score_avg"])) == expected[key] / 4
    with gzip.open(DATA / "pgsc_calc_v2.3.0/aggregated_scores.txt.gz", "rt") as handle:
        aggregate = list(csv.DictReader(handle, delimiter="\t"))
    assert len(aggregate) == 6
    assert {(r["PGS"], r["IID"]) for r in aggregate} == set(expected)
    for row in aggregate:
        key = row["PGS"], row["IID"]
        assert Decimal(row["SUM"]) == reference["1", *key]["score"]
        assert Decimal(row["DENOM"]) == reference["1", *key]["denom"]
        assert Decimal(row["AVG"]) == expected[key] / 4
