"""Cross-score routing and hard-call arithmetic, independently checked with Decimal."""

import csv
import gzip
import runpy
from collections import defaultdict
from decimal import Decimal
from pathlib import Path

import polars as pl
import pysam

from .parity_helpers import read_reference, run_case, validate_inputs

CASE = Path(__file__).parent.parent / "data/v2_parity/arithmetic"
SCORES = sorted(CASE.glob("PGS*.txt.gz"))
TARGETS = [(CASE / f"target.chr{chrom}.vcf.gz", chrom) for chrom in ("1", "2")]


def exact_scores():
    """Count alleles from GT directly, without using v3 matching/dosage code."""
    calls = {}
    for path, chrom in TARGETS:
        with pysam.VariantFile(path) as vcf:
            for variant in vcf:
                for sample, call in variant.samples.items():
                    calls[chrom, variant.pos, sample] = tuple(
                        variant.alleles[index] for index in call["GT"]
                    )
    expected = defaultdict(lambda: defaultdict(Decimal))
    identities = {}
    for path in SCORES:
        with gzip.open(path, "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                chrom, position = row["chr_name"], int(row["chr_position"])
                identities[row["accession"], int(row["row_nr"])] = (
                    chrom,
                    position,
                    row["effect_allele"],
                )
                for sample in ("S0", "S1", "S2", "S3"):
                    gt = calls[chrom, position, sample]
                    dosage = gt.count(row["effect_allele"])
                    if row["effect_type"] == "dominant":
                        dosage = int(dosage > 0)
                    elif row["effect_type"] == "recessive":
                        dosage = int(dosage == 2)
                    values = expected[chrom, row["accession"], sample]
                    values["score"] += dosage * Decimal(row["effect_weight"])
                    values["dosage_sum"] += dosage
                    values["allele_count"] += len(gt)
    return expected, identities


def reference_budgets():
    """Accumulate six-significant-figure rounding bounds before cancellation."""
    budgets = defaultdict(Decimal)
    for path in (CASE / "pgsc_calc_v2.3.0").glob("*.sscore"):
        chrom = path.stem.rsplit("_", 3)[1]
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                for column, text in row.items():
                    if column.startswith("PGS") and column.endswith("_SUM"):
                        value = Decimal(text)
                        # Exact zero in these fully called files needs no text
                        # rounding allowance. Nonzero shards may cancel later.
                        if value:
                            budgets[chrom, column.removesuffix("_SUM"), row["IID"]] += (
                                Decimal("0.5") * Decimal(10) ** (value.adjusted() - 5)
                            )
    return budgets


def test_arithmetic_inputs_regenerate(tmp_path):
    """The generator cannot silently change the committed reference inputs."""
    runpy.run_path(str(CASE / "generate.py"))["generate"](tmp_path)
    generated = list(tmp_path.iterdir())
    assert len(generated) == 12
    for path in generated:
        assert path.read_bytes() == (CASE / path.name).read_bytes(), path.name


def test_multiple_scores_decimal_weights_and_effect_models(tmp_path, monkeypatch):
    validate_inputs(CASE)
    expected, identities = exact_scores()
    reference = read_reference(CASE / "pgsc_calc_v2.3.0")
    assert set(reference) == set(expected)
    # PLINK's text output has six significant figures. One shard can therefore
    # differ by at most half a unit at that precision; the public v3 CSV adds
    # half a unit at its six decimal places. The bound is defined from precision,
    # never fitted to differences observed between the two implementations.
    reference_budget = reference_budgets()
    for key, values in reference.items():
        assert abs(values["score"] - expected[key]["score"]) <= reference_budget[key], (
            key
        )
        assert values["denom"] == values["allele_count"]
        # v2 batches scores together and reports union-of-variant statistics in
        # each .sscore. These are not accession-specific denominators/dosages.
        # PGS990004 has one contributing variant per effect-type shard, so its
        # raw statistics are comparable; the others intentionally are not.
        if key[1] == "PGS990004":
            for field in ("allele_count", "dosage_sum"):
                assert values[field] == expected[key][field], (key, field)
        else:
            # additive_0 names A,C,A,C,C,C at variant indices 1..6: counting
            # those GT alleles gives (5,9,4,3) in S0..S3 order on chr1.
            # PGS990001 also appears in additive_1, which names C at position
            # 1000 and adds (0,1,2,1). Chr2 reverses the sample call vectors.
            # Pin raw union statistics as well as the per-score v3 statistics;
            # merely asserting that the denominators differ misses corruption.
            dosage_counts = (5, 10, 6, 4) if key[1] == "PGS990001" else (5, 9, 4, 3)
            if key[0] == "2":
                dosage_counts = dosage_counts[::-1]
            assert (
                values["dosage_sum"]
                == dict(zip(("S0", "S1", "S2", "S3"), dosage_counts, strict=True))[
                    key[2]
                ]
            ), key
            assert values["allele_count"] == (14 if key[1] == "PGS990001" else 12), key

    # Six cached variants and batch_size=2 exercise real Dask/Zarr boundaries
    # while keeping the fixture inspectable. Opposing weights straddle chunks
    # in the expanded effect-allele/model matrix, where routing errors matter.
    actual = run_case(TARGETS, SCORES, tmp_path, monkeypatch, batch_size=2)
    for index, (_, chrom) in enumerate(TARGETS):
        rows = pl.read_parquet(
            tmp_path / "score" / str(index) / "partial_scores.parquet"
        ).to_dicts()
        keys = {(chrom, row["accession"], row["sample_id"]) for row in rows}
        assert len(rows) == len(keys) == 16
        assert keys == {key for key in expected if key[0] == chrom}
        for row in rows:
            key = chrom, row["accession"], row["sample_id"]
            for field in ("score", "allele_count", "dosage_sum"):
                value = Decimal(str(row[field]))
                # Decimal expectations also catch a shared error in the v2/v3
                # comparison. Float64 arithmetic permits tiny round-off only.
                budget = Decimal("0.000000001") if field == "score" else 0
                assert abs(value - expected[key][field]) <= budget, (key, field)
                if field == "score":
                    assert abs(value - reference[key][field]) <= reference_budget[
                        key
                    ] + Decimal("0.000000001"), key

    totals = defaultdict(lambda: defaultdict(Decimal))
    for (_, accession, sample), values in expected.items():
        for field, value in values.items():
            totals[accession, sample][field] += value
    rows = actual.to_dicts()
    assert len(rows) == len(totals) == 16
    assert {(row["accession"], row["sample_id"]) for row in rows} == set(totals)
    for row in rows:
        key = row["accession"], row["sample_id"]
        for field in ("score", "allele_count", "dosage_sum", "score_avg"):
            expected_value = (
                totals[key]["score"] / totals[key]["allele_count"]
                if field == "score_avg"
                else totals[key][field]
            )
            budget = Decimal("0.000000500001") if field in ("score", "score_avg") else 0
            assert abs(Decimal(str(row[field])) - expected_value) <= budget, (
                key,
                field,
            )
        assert int(row["n_matched"]) == sum(
            accession == key[0] for accession, _ in identities
        )

    # v2's aggregator sums every shard's denominator for every accession,
    # including shards with no column for that score. Preserve that evidence,
    # but assert v3's per-score biological denominator against the GT oracle.
    with gzip.open(CASE / "pgsc_calc_v2.3.0/aggregated_scores.txt.gz", "rt") as handle:
        aggregate = list(csv.DictReader(handle, delimiter="\t"))
    assert len(aggregate) == len(totals)
    assert {(r["PGS"], r["IID"]) for r in aggregate} == set(totals)
    for row in aggregate:
        key = row["PGS"], row["IID"]
        saved_sum = sum(reference[chrom, *key]["score"] for chrom in ("1", "2"))
        assert abs(Decimal(row["SUM"]) - saved_sum) < Decimal("0.000000001")
        assert Decimal(row["DENOM"]) == 36
        assert abs(Decimal(row["AVG"]) - saved_sum / 36) < Decimal("0.000000001")
        assert Decimal(row["DENOM"]) != totals[key]["allele_count"]

    with gzip.open(CASE / "pgsc_calc_v2.3.0/synthetic_log.csv.gz", "rt") as handle:
        selected_v2 = [
            row for row in csv.DictReader(handle) if row["match_status"] == "matched"
        ]
    selected_v3 = []
    for path in (tmp_path / "matched/logs").rglob("*.csv.gz"):
        with gzip.open(path, "rt") as handle:
            selected_v3.extend(
                row
                for row in csv.DictReader(handle)
                if row["match_summary"] == "matched"
            )
    for selected in (selected_v2, selected_v3):
        observed = {
            (row["accession"], int(row["row_nr"])): (
                row["chr_name"],
                int(row["chr_position"]),
                row["matched_effect_allele"],
            )
            for row in selected
        }
        assert len(selected) == len(observed) == len(identities)
        assert observed == identities
