"""Score identities survive input order and complete hard-call sample batches."""

import gzip
from collections import defaultdict
from decimal import Decimal

import polars as pl
import pysam
import pytest

from .parity_helpers import read_reference, run_case, validate_inputs
from .test_v2_parity_arithmetic import (
    CASE,
    SCORES,
    TARGETS,
    exact_scores,
    reference_budgets,
)


@pytest.mark.parametrize("layout", ["sample_order", "input_order", "sample_batches"])
def test_arithmetic_partition_invariance(layout, tmp_path, monkeypatch):
    validate_inputs(CASE)
    targets = list(TARGETS)
    scores = list(SCORES)
    if layout == "input_order":
        targets.reverse()
        scores.reverse()
    else:
        targets = []
        for source, chrom in TARGETS:
            # Permute columns explicitly: selecting samples through a VCF reader
            # may preserve header order. These transformations retain the exact
            # called GT values and do not change an imputation population.
            with gzip.open(source, "rt") as handle:
                lines = handle.readlines()
            header = next(line for line in lines if line.startswith("#CHROM"))
            samples = header.rstrip().split("\t")[9:]
            groups = (
                [("S2", "S0"), ("S3", "S1")]
                if layout == "sample_batches"
                else [
                    ("S3", "S1", "S0", "S2")
                    if chrom == "1"
                    else ("S1", "S2", "S3", "S0")
                ]
            )
            for index, group in enumerate(groups):
                columns = [*range(9), *(9 + samples.index(sample) for sample in group)]
                plain = tmp_path / f"chr{chrom}.{index}.vcf"
                plain.write_text(
                    "".join(
                        line
                        if line.startswith("##")
                        else "\t".join(
                            line.rstrip().split("\t")[column] for column in columns
                        )
                        + "\n"
                        for line in lines
                    )
                )
                compressed = plain.with_suffix(".vcf.gz")
                pysam.tabix_compress(str(plain), str(compressed), force=True)
                pysam.tabix_index(str(compressed), preset="vcf", force=True)
                targets.append((compressed, chrom))

    actual = run_case(targets, scores, tmp_path / "run", monkeypatch, batch_size=2)
    exact, _ = exact_scores()
    reference = read_reference(CASE / "pgsc_calc_v2.3.0")
    assert set(reference) == set(exact)
    budgets = reference_budgets()
    seen = set()
    for index, (_, chrom) in enumerate(targets):
        rows = pl.read_parquet(
            tmp_path / "run/score" / str(index) / "partial_scores.parquet"
        ).to_dicts()
        for row in rows:
            key = chrom, row["accession"], row["sample_id"]
            assert key not in seen
            seen.add(key)
            for field in ("score", "allele_count", "dosage_sum"):
                tolerance = Decimal("0.000000001") if field == "score" else 0
                assert abs(Decimal(str(row[field])) - exact[key][field]) <= tolerance, (
                    layout,
                    key,
                    field,
                )
            # Saved text rounds to six significant figures; use the same
            # shardwise precision bound as the arithmetic comparison.
            assert abs(Decimal(str(row["score"])) - reference[key]["score"]) <= budgets[
                key
            ] + Decimal("0.000000001"), key
    assert seen == set(exact)
    totals = defaultdict(lambda: defaultdict(Decimal))
    for (_, accession, sample), values in exact.items():
        for field, value in values.items():
            totals[accession, sample][field] += value
    rows = actual.to_dicts()
    assert len(rows) == len(totals) == 16
    assert {(r["accession"], r["sample_id"]) for r in rows} == set(totals)
    for row in rows:
        values = totals[row["accession"], row["sample_id"]]
        for field in ("score", "allele_count", "dosage_sum", "score_avg"):
            expected = (
                values["score"] / values["allele_count"]
                if field == "score_avg"
                else values[field]
            )
            tolerance = (
                Decimal("0.000000500001") if field in ("score", "score_avg") else 0
            )
            assert abs(Decimal(str(row[field])) - expected) <= tolerance, (
                layout,
                row,
                field,
            )
