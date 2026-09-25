"""Match policy affects biological allele identity before any weighted arithmetic.

Compare shared policies with the unchanged v2 workflow; pin deliberate v3 changes
with independently counted GT alleles rather than preserving old matcher defects.
"""

import csv
import gzip
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import polars as pl
import pysam
import pytest

from .parity_helpers import read_reference, run_case, validate_inputs

CASE = Path(__file__).parent.parent / "data/v2_parity/matching"
MODES = {
    "default": (),
    "keep_ambiguous": ("--keep_ambiguous",),
    "ignore_strand_flips": ("--ignore_strand_flips",),
}
# The expected target record and effect sequence are explicit biological choices,
# not inferred by calling the matcher or reading its output.
SELECTED = {
    "PGS990020": {
        0: ("direct", "A"),
        1: ("swapped", "A"),
        2: ("flip_ref", "T"),
        3: ("flip_alt", "T"),
        4: ("no_other", "C"),
    },
    "PGS990021": {0: ("anchor", "C")},
    "PGS990022": {
        0: ("anchor", "C"),
        1: ("distinct_alt", "C"),
        2: ("distinct_alt", "G"),
    },
    "PGS990023": {0: ("anchor", "C"), 1: ("candidate_ref", "A")},
    "PGS990024": {
        0: ("anchor", "C"),
        1: ("reverse_complement", "GT"),
        4: ("forward_indel", "AT"),
    },
    "PGS990025": {0: ("anchor", "C")},
}


def expected_selection(mode):
    selected = {accession: dict(rows) for accession, rows in SELECTED.items()}
    if mode == "keep_ambiguous":
        selected["PGS990025"].update(
            {1: ("palindrome", "A"), 2: ("no_other_ambiguous", "C")}
        )
    if mode == "ignore_strand_flips":
        for row in (2, 3):
            del selected["PGS990020"][row]
        del selected["PGS990024"][1]
    return selected


def hardcall_expectations(selected):
    with pysam.VariantFile(CASE / "target.chr1.vcf.gz") as vcf:
        samples = list(vcf.header.samples)
        records = {record.id: record for record in vcf}
    expected = {}
    for accession, choices in selected.items():
        with gzip.open(CASE / f"{accession}.txt.gz", "rt") as handle:
            weights = {
                int(r["row_nr"]): Decimal(r["effect_weight"])
                for r in csv.DictReader(handle, delimiter="\t")
            }
        for sample in samples:
            total, dosage = Decimal(0), 0
            for row, (record_id, effect) in choices.items():
                record = records[record_id]
                copies = record.samples[sample]["GT"].count(
                    record.alleles.index(effect)
                )
                total += copies * weights[row]
                dosage += copies
            # Each retained score row observes two alleles. Distinct ALT effects
            # at one locus are separate score rows, not a collision (ADR 0010).
            count = 2 * len(choices)
            expected[accession, sample] = {
                "score": total,
                "dosage_sum": Decimal(dosage),
                "allele_count": Decimal(count),
                "score_avg": (total / count).quantize(
                    Decimal("0.000001"), rounding=ROUND_HALF_UP
                ),
            }
    return expected


@pytest.mark.parametrize("mode", MODES)
def test_matching_policy_and_scores(mode, tmp_path, monkeypatch):
    selected = expected_selection(mode)
    expected = hardcall_expectations(selected)
    result = run_case(
        [(CASE / "target.chr1.vcf.gz", "1")],
        sorted(CASE.glob("PGS*.txt.gz")),
        tmp_path,
        monkeypatch,
        match_flags=("--keep_multiallelic", *MODES[mode]),
    )
    actual = result.to_dicts()
    assert len(actual) == len(expected) == 18
    assert {(r["accession"], r["sample_id"]) for r in actual} == set(expected)
    for row in actual:
        key = row["accession"], row["sample_id"]
        assert row["n_matched"] == len(selected[row["accession"]])
        for field, value in expected[key].items():
            assert Decimal(str(row[field])) == value, (mode, key, field)

    partial = pl.read_parquet(tmp_path / "score/0/partial_scores.parquet").to_dicts()
    assert len(partial) == len(expected)
    assert {(r["accession"], r["sample_id"]) for r in partial} == set(expected)
    for row in partial:
        key = row["accession"], row["sample_id"]
        for field in ("score", "dosage_sum", "allele_count"):
            assert Decimal(str(row[field])) == expected[key][field], (key, field)

    # Validate actual genotype routing, not just the allele string reconstructed
    # in logs. This VCF's only absent score position lies after every real row.
    with pysam.VariantFile(CASE / "target.chr1.vcf.gz") as vcf:
        records = list(vcf)
    table = pl.read_parquet(tmp_path / "matched/match_table.parquet").to_dicts()
    actual_selected = {}
    for row in table:
        if not row["match_result"]["is_matched"]:
            continue
        record = records[row["target_row_nr"]]
        effect = record.alleles[row["match_result"]["effect_allele_idx"]]
        key = row["accession"], row["row_nr"]
        assert key not in actual_selected
        actual_selected[key] = record.id, effect
    assert actual_selected == {
        (a, n): identity for a, rows in selected.items() for n, identity in rows.items()
    }

    logs = []
    for path in (tmp_path / "matched/logs").rglob("*.csv.gz"):
        with gzip.open(path, "rt") as handle:
            logs.extend(csv.DictReader(handle))
    outcomes = {
        (r["accession"], int(r["row_nr"])): (r["match_summary"], r["exclusion_reason"])
        for r in logs
    }
    assert len(logs) == len(outcomes) == 23
    for row in (5, 6):
        assert outcomes["PGS990020", row][0] == "unmatched"
    for row in (1, 2):
        assert outcomes["PGS990021", row] == ("excluded", "collision")
    for row in (2, 3):
        assert outcomes["PGS990024", row][0] == "unmatched"
    if mode != "keep_ambiguous":
        for row in (1, 2):
            assert outcomes["PGS990025", row] == ("excluded", "ambiguous")
    if mode == "ignore_strand_flips":
        for accession, row in (("PGS990020", 2), ("PGS990020", 3), ("PGS990024", 1)):
            assert outcomes[accession, row] == ("excluded", "strand_flip")

    reference = CASE / "pgsc_calc_v2.3.0" / mode
    validate_inputs(CASE, reference)
    saved = read_reference(reference)
    # Pin the old policy with explicit selections, including its deliberate
    # differences: ID-based collisions, plain complements and indel inference.
    v2_selected = expected_selection(mode)
    v2_selected["PGS990022"] = {0: ("anchor", "C")}
    v2_selected["PGS990024"] = {
        0: ("anchor", "C"),
        2: ("plain_complement", "TG"),
        3: ("indel_inference", "T"),
        4: ("forward_indel", "AT"),
    }
    if mode == "ignore_strand_flips":
        del v2_selected["PGS990024"][2]
        del v2_selected["PGS990024"][3]
    v2_selected["PGS990025"][2] = ("no_other_ambiguous", "C")
    v2_expected = hardcall_expectations(v2_selected)
    by_id = {record.id: record for record in records}
    with gzip.open(reference / "synthetic_log.csv.gz", "rt") as handle:
        chosen = [r for r in csv.DictReader(handle) if r["match_status"] == "matched"]
    identities = {
        (r["accession"], int(r["row_nr"])): (
            int(r["chr_position"]),
            r["REF"],
            r["matched_effect_allele"],
        )
        for r in chosen
    }
    wanted_identities = {
        (a, n): (by_id[vid].pos, by_id[vid].ref, effect)
        for a, rows in v2_selected.items()
        for n, (vid, effect) in rows.items()
    }
    assert len(chosen) == len(identities)
    assert identities == wanted_identities
    assert set(saved) == {("1", *key) for key in v2_expected}
    # PLINK scores these accessions together: its count/dosage columns describe
    # the UNION of selected alleles, including rows with zero-filled weights in
    # other score columns. v3 reports per-score statistics. This is an explicit
    # denominator difference, not a tolerance or a reason to skip AVG checks.
    union = {identity for rows in v2_selected.values() for identity in rows.values()}
    for key, values in v2_expected.items():
        sample = key[1]
        observed = saved["1", *key]
        assert observed["score"] == values["score"], (mode, key)
        assert observed["allele_count"] == observed["denom"] == 2 * len(union)
        assert observed["dosage_sum"] == sum(
            by_id[vid].samples[sample]["GT"].count(by_id[vid].alleles.index(effect))
            for vid, effect in union
        )
    with gzip.open(reference / "aggregated_scores.txt.gz", "rt") as handle:
        aggregate = list(csv.DictReader(handle, delimiter="\t"))
    assert len(aggregate) == len(v2_expected)
    assert {(r["PGS"], r["IID"]) for r in aggregate} == set(v2_expected)
    for row in aggregate:
        key = row["PGS"], row["IID"]
        raw = saved["1", *key]
        assert Decimal(row["SUM"]) == raw["score"]
        assert Decimal(row["DENOM"]) == raw["denom"]
        # v2's aggregate is written to more precision than v3's six decimals.
        quantum = Decimal("0.000001")
        assert Decimal(row["AVG"]).quantize(quantum, rounding=ROUND_HALF_UP) == (
            raw["score"] / raw["denom"]
        ).quantize(quantum, rounding=ROUND_HALF_UP)
    # These three accessions share v2's rules, including exclusion of colliding
    # score rows and priority selection between multiple target candidates.
    for accession in ("PGS990020", "PGS990021", "PGS990023"):
        for sample in ("S0", "S1", "S2"):
            key = accession, sample
            assert saved["1", *key]["score"] == expected[key]["score"], (mode, key)


def test_matching_overlap_rejects_incomplete_accessions(tmp_path, monkeypatch):
    """A high overlap gate must remove scores, not report their partial sums.

    Only the distinct-ALT and candidate-selection scores match every row. The
    former deliberately differs from v2's ID-based collision policy (ADR 0010).
    """
    result = run_case(
        [(CASE / "target.chr1.vcf.gz", "1")],
        sorted(CASE.glob("PGS*.txt.gz")),
        tmp_path,
        monkeypatch,
        match_flags=("--keep_multiallelic",),
        min_overlap="0.8",
    )
    wanted = {"PGS990022", "PGS990023"}
    selected = {a: rows for a, rows in SELECTED.items() if a in wanted}
    expected = hardcall_expectations(selected)
    assert result.height == 6
    assert {(r["accession"], r["sample_id"]) for r in result.to_dicts()} == set(
        expected
    )
    for row in result.to_dicts():
        key = row["accession"], row["sample_id"]
        for field, value in expected[key].items():
            assert Decimal(str(row[field])) == value, (key, field)


def test_matching_inputs_are_reproducible(tmp_path):
    import runpy

    runpy.run_path(str(CASE / "generate.py"))["generate"](tmp_path)
    generated = sorted(p for p in tmp_path.iterdir() if p.is_file())
    assert len(generated) == 14
    for path in generated:
        assert path.read_bytes() == (CASE / path.name).read_bytes(), path.name
