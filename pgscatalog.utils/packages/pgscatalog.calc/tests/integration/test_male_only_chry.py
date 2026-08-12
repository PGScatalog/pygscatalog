"""Acceptance tests for chrY from a male-only file (ADR 0007).

Y is the only chromosome a person can lack entirely, so a chrY leaf covering a
subset of the sampleset is biology rather than a truncated file. That relaxes one
of ADR 0003's four cache shape invariants - cross-chromosome equality of the sample
union - for Y and for nothing else.

1000G ships no chrY genotypes, so the fixture is real non-pseudoautosomal chrX data
with the women dropped and the contig renamed. What these tests need from it is its
shape: hemizygous calls, and a sample list that is a strict subset of the
sampleset's.
"""

import gzip
import pathlib

import pysam
import pytest

from .conftest import (
    _is_ambiguous,
    _subset_samples,
    _write_vcf,
    load_targets,
    run_cli,
    run_serial,
)

# enough chrY variants to be worth a sum, few enough to keep the fixture quick
N_Y_VARIANTS = 10


def _observed(records: list) -> tuple[dict[str, int], dict[str, int]]:
    """What the file says: each sample's ALT dosage and allele count"""
    dosage: dict[str, int] = {}
    alleles: dict[str, int] = {}
    for rec in records:
        for name, call in rec.samples.items():
            gt = call["GT"]
            dosage[name] = dosage.get(name, 0) + sum(allele == 1 for allele in gt)
            alleles[name] = alleles.get(name, 0) + len(gt)
    return dosage, alleles


@pytest.fixture(scope="session")
def chry_target(vcf_path, tmp_path_factory) -> pathlib.Path:
    """A male-only chrY file: every call in it is hemizygous, so nothing about it
    is ambiguous except which samples are absent"""
    out_dir = tmp_path_factory.mktemp("chry")

    males: set[str] = set()
    records: list = []
    with pysam.VariantFile(str(vcf_path)) as vcf:
        header = vcf.header.copy()
        for rec in vcf.fetch("X"):
            if len(rec.ref) != 1 or len(rec.alts) != 1 or len(rec.alts[0]) != 1:
                continue
            if _is_ambiguous(rec.ref, rec.alts[0]):
                continue
            hemizygous = {s for s, call in rec.samples.items() if len(call["GT"]) == 1}
            if not hemizygous:
                continue  # pseudoautosomal, so nobody here is hemizygous
            males = males or hemizygous
            if hemizygous == males:
                records.append(rec)
            if len(records) == N_Y_VARIANTS:
                break

    assert len(records) == N_Y_VARIANTS, "not enough hemizygous test variants"

    x_path = out_dir / "chrX.vcf.gz"
    _write_vcf(x_path, header, records)
    male_only = _subset_samples(x_path, out_dir / "male_only.vcf.gz", sorted(males))

    rename = out_dir / "rename.txt"
    rename.write_text("X\tY\n")
    path = out_dir / "chrY.vcf.gz"
    pysam.bcftools.annotate(
        "--rename-chrs",
        str(rename),
        "-Oz",
        "-o",
        str(path),
        str(male_only),
        catch_stdout=False,
    )
    pysam.bcftools.index(str(path))
    return path


@pytest.fixture(scope="session")
def chry_records(chry_target) -> list:
    """The male-only file as it was written, which is what the pipeline reads"""
    with pysam.VariantFile(str(chry_target)) as vcf:
        return list(vcf.fetch("Y"))


@pytest.fixture(scope="session")
def autosome_and_chry_scorefile(
    target_records, chry_records, tmp_path_factory
) -> pathlib.Path:
    """One accession spanning chromosome 1 and chrY, weight 1 per variant so the
    score is the effect allele dosage sum and the acceptance numbers are readable"""
    _, by_chrom = target_records
    path = tmp_path_factory.mktemp("chry_score") / "PGS999993.txt.gz"

    rows = [("1", rec) for rec in by_chrom["1"]] + [("Y", r) for r in chry_records]
    with gzip.open(path, "wt") as f:
        f.write(
            "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\t"
            "effect_type\tis_duplicated\taccession\trow_nr\n"
        )
        for row_nr, (chrom, rec) in enumerate(rows):
            f.write(
                f"{chrom}\t{rec.pos}\t{rec.alts[0]}\t{rec.ref}\t1.0\tadditive\t"
                f"False\tPGS999993\t{row_nr}\n"
            )

    return path


def test_a_male_only_chry_passes_match_combine_and_scores(
    chromosome_targets,
    target_records,
    chry_target,
    chry_records,
    autosome_and_chry_scorefile,
    tmp_path,
    monkeypatch,
):
    """Autosomal leaves cover everyone, the chrY leaf covers the men, and the run
    completes. Women get no chrY contribution and an allele count reduced by
    exactly the matched chrY variants times their ploidy there - which is what
    SCORE_COMBINE's `GROUP BY (sampleset, accession, sample_id)` and `SUM` already
    produce for a sample no chrY leaf holds (ADR 0003, corrected by ADR 0007).
    """
    _, by_chrom = target_records
    chr1_dosage, chr1_alleles = _observed(by_chrom["1"])
    y_dosage, y_alleles = _observed(chry_records)

    scores = run_serial(
        [chromosome_targets[0], (chry_target, "Y")],
        tmp_path,
        autosome_and_chry_scorefile,
        monkeypatch,
    )

    # everyone gets exactly one row, including the women no chrY leaf holds
    assert scores.shape[0] == 3202
    assert scores["sample_id"].n_unique() == 3202
    # n_matched is genome-wide and sample-independent, so chrY counts for everyone
    assert scores["n_matched"].unique().to_list() == [len(by_chrom["1"]) + N_Y_VARIANTS]

    for row in scores.iter_rows(named=True):
        sample = row["sample_id"]
        assert row["dosage_sum"] == chr1_dosage[sample] + y_dosage.get(sample, 0)
        assert row["allele_count"] == chr1_alleles[sample] + y_alleles.get(sample, 0)

    men = set(y_alleles)
    women = set(chr1_alleles) - men
    assert len(men) == 1599 and len(women) == 1603

    # men are hemizygous on Y: one allele per matched chrY variant, never two
    assert set(y_alleles.values()) == {N_Y_VARIANTS}
    assert max(y_dosage.values()) <= N_Y_VARIANTS

    # a woman's allele count is short by the matched chrY variants times the ploidy
    # a sample covered there carries, which is 1
    counts = dict(zip(scores["sample_id"], scores["allele_count"], strict=True))
    assert {counts[s] for s in women} == {2 * len(by_chrom["1"])}
    assert {counts[s] for s in men} == {2 * len(by_chrom["1"]) + N_Y_VARIANTS * 1}

    # and no chrY contribution at all, so her score is her autosomal one
    dosages = dict(zip(scores["sample_id"], scores["dosage_sum"], strict=True))
    assert all(dosages[s] == chr1_dosage[s] for s in women)
    assert scores["score"].to_list() == pytest.approx(scores["dosage_sum"].to_list())


def test_a_truncated_chry_is_accepted_like_a_male_only_one(
    chromosome_targets,
    target_records,
    chry_target,
    chry_records,
    autosome_and_chry_scorefile,
    tmp_path,
    monkeypatch,
):
    """The known cost of the carve-out, asserted deliberately so nobody later reads
    it as a bug: a chrY file that lost samples is byte-for-byte the shape of one
    that never had them, so the pipeline cannot refuse it"""
    _, by_chrom = target_records
    with pysam.VariantFile(str(chry_target)) as vcf:
        kept = sorted(vcf.header.samples)[:100]

    truncated = _subset_samples(chry_target, tmp_path / "truncated_chrY.vcf.gz", kept)

    scores = run_serial(
        [chromosome_targets[0], (truncated, "Y")],
        tmp_path,
        autosome_and_chry_scorefile,
        monkeypatch,
    )

    assert scores.shape[0] == 3202
    counts = dict(zip(scores["sample_id"], scores["allele_count"], strict=True))
    covered = {s for s, n in counts.items() if n > 2 * len(by_chrom["1"])}
    assert covered == set(kept)


def test_a_truncated_autosome_is_still_rejected(
    chromosome_targets, scorefile, tmp_path, monkeypatch
):
    """The relaxation must not leak: every chromosome except Y still has to cover
    the sampleset's full sample union exactly"""
    full_chr8, chrom = chromosome_targets[1]
    with pysam.VariantFile(str(full_chr8)) as vcf:
        kept = sorted(vcf.header.samples)[:100]

    truncated = _subset_samples(full_chr8, tmp_path / "truncated_chr8.vcf.gz", kept)
    targets = [chromosome_targets[0], (truncated, chrom)]

    zips = load_targets(targets, tmp_path / "load", scorefile, monkeypatch)
    samples = []
    for (_, c), zip_path in zip(targets, zips, strict=True):
        out_dir = tmp_path / "match" / c
        run_cli(
            [
                "pgsc_calc",
                "match",
                "--zarr_zip_file",
                str(zip_path),
                "--score_paths",
                str(scorefile),
                "--out_dir",
                str(out_dir),
            ],
            monkeypatch,
        )
        samples.append(out_dir / "samples.parquet")

    with pytest.raises(ValueError) as excinfo:
        run_cli(
            [
                "pgsc_calc",
                "match-combine",
                "--matches",
                *[str(p.parent / "matches.parquet") for p in samples],
                "--samples",
                *[str(p) for p in samples],
                "--score_paths",
                str(scorefile),
                "--sampleset",
                "test",
                "--out_dir",
                str(tmp_path / "match_combine"),
            ],
            monkeypatch,
        )

    message = str(excinfo.value)
    assert f"['{chrom}']" in message
    assert "must cover the same samples" in message


def test_a_mixed_sex_chry_vcf_is_rejected(
    chry_target, chry_records, autosome_and_chry_scorefile, tmp_path, monkeypatch
):
    """VCF writes `.` for a no-call and for a sample that has no Y alike, so a
    mixed-sex chrY file is uninterpretable rather than imputable (ADR 0007).
    The message names the format as the limit; it must not suggest reshaping the
    data into a male-only file, which would be a wrong answer for anyone whose
    file genuinely is mixed-sex (ADR 0004's precedent).
    """
    with pysam.VariantFile(str(chry_target)) as vcf:
        header = vcf.header.copy()
        records = list(vcf.fetch("Y"))

    # `.`, spelled as pysam wants it: a call of ploidy one with no allele in it.
    # In a male-only file that is a hemizygous no-call, but nothing in a chrY file
    # says whether the sample has a Y at all
    for call in list(records[0].samples.values())[:10]:
        call["GT"] = (None,)

    path = tmp_path / "mixed_sex_chrY.vcf.gz"
    _write_vcf(path, header, records)

    with pytest.raises(ValueError) as excinfo:
        load_targets(
            [(path, "Y")], tmp_path / "cache", autosome_and_chry_scorefile, monkeypatch
        )

    message = str(excinfo.value)
    assert "VCF" in message
    assert "male-only" not in message.lower()
