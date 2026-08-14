"""Acceptance test for haploid calls from bgen, end to end (ADR 0007).

The VCF path is the reference implementation, so the thing worth asserting is that
the same chrX calls score the same whichever format they arrive in. Ploidy reaches
the bgen reader on a different channel (``BgenVar.ploidy``) than it reaches the VCF
reader (the length of ``GT``), and the cache is supposed to be blind to that.
"""

import gzip
import pathlib
import subprocess

import numpy as np
import polars as pl
import pysam
import pytest
from bgen import BgenWriter
from polars.testing import assert_frame_equal

from .conftest import SAMPLESET, _is_ambiguous, _write_vcf, load_target, run_cli


@pytest.fixture(scope="session")
def big_bgen_n_samples():
    return 10_000


@pytest.fixture(scope="session")
def big_bgen_n_variants():
    return 500_000


@pytest.fixture(scope="session")
def really_big_bgen(tmp_path_factory, big_bgen_n_samples, big_bgen_n_variants):
    base_dir = tmp_path_factory.mktemp("big_file", numbered=False)
    bgen_path = base_dir / "test_large_file.bgen"
    bgi_path = base_dir / "test_large_file.bgen.bgi"

    # this file ~1GB
    rng = np.random.default_rng()
    with BgenWriter(bgen_path, n_samples=big_bgen_n_samples) as bfile:
        for i in range(big_bgen_n_variants):
            geno = rng.random((big_bgen_n_samples, 3)).astype(np.float16)
            bfile.add_variant(
                varid=f"var{i}",
                rsid="rsid{i}",
                chrom=f"chr{i}",
                pos=i,
                alleles=["A", "G"],
                genotypes=geno,
                bit_depth=1,
            )

    sample_ids = [f"sample{i}" for i in range(big_bgen_n_samples)]

    return sample_ids, bgen_path, bgi_path


def _probabilities(rec, n_samples):
    """Write each hard call as a probability of 1, so the bgen holds exactly the
    calls the VCF does. A haploid row is two wide and a diploid one three, padded
    with NaN in the trailing slot the sample doesn't have
    """
    ploidy = np.array([len(x["GT"]) for x in rec.samples.values()], dtype=np.uint8)
    probs = np.full((n_samples, int(ploidy.max()) + 1), np.nan)
    for i, call in enumerate(rec.samples.values()):
        gt = call["GT"]
        probs[i, : len(gt) + 1] = 0.0
        probs[i, sum(gt)] = 1.0  # unphased: the column is the ALT allele count
    return probs, ploidy


@pytest.fixture(scope="session")
def chrx_records(vcf_path):
    """Fully called biallelic chrX SNPs: 1000G writes a man's non-pseudoautosomal
    call as `1` and a woman's as `0/1`, so ploidy is mixed within every variant.

    Missing calls are dropped, because VCF can't say whether `.` is a haploid
    no-call or a diploid one while bgen can - a difference the two readers are
    meant to have, and one this comparison isn't about
    """
    with pysam.VariantFile(str(vcf_path)) as vcf:
        header = vcf.header.copy()
        records = [
            rec
            for rec in vcf.fetch("X")
            if len(rec.ref) == 1
            and len(rec.alts) == 1
            and len(rec.alts[0]) == 1
            and not _is_ambiguous(rec.ref, rec.alts[0])
            and all(None not in call["GT"] for call in rec.samples.values())
        ][:20]

    assert len(records) > 10, "not enough chrX test variants"
    return header, records


@pytest.fixture(scope="session")
def chrx_vcf_target(chrx_records, tmp_path_factory) -> pathlib.Path:
    header, records = chrx_records
    path = tmp_path_factory.mktemp("chrx_vcf") / "chrX.vcf.gz"
    _write_vcf(path, header, records)
    return path


@pytest.fixture(scope="session")
def chrx_bgen_target(chrx_records, tmp_path_factory) -> tuple[pathlib.Path, ...]:
    header, records = chrx_records
    out_dir = tmp_path_factory.mktemp("chrx_bgen")
    path = out_dir / "chrX.bgen"
    samples = list(header.samples)

    with BgenWriter(path, n_samples=len(samples), samples=samples) as bfile:
        for i, rec in enumerate(records):
            probs, ploidy = _probabilities(rec, len(samples))
            bfile.add_variant(
                varid=f"var{i}",
                rsid=f"rs{i}",
                chrom="X",
                pos=rec.pos,
                alleles=[rec.ref, rec.alts[0]],
                genotypes=probs,
                ploidy=ploidy,
                phased=False,
                bit_depth=8,
            )

    # BgenWriter leaves its own .bgi behind, but bgenix won't query one it didn't
    # write, so overwrite it
    subprocess.run(["bgenix", "-index", "-clobber", "-g", str(path)], check=True)
    sample_path = out_dir / "chrX.sample"
    sample_path.write_text(
        "ID_1 ID_2 missing\n0 0 0\n" + "".join(f"{x} {x} 0\n" for x in samples)
    )
    return path, sample_path


@pytest.fixture(scope="session")
def chry_records(chrx_records) -> tuple[list[str], list]:
    """The non-pseudoautosomal chrX variants, where the same men are hemizygous
    everywhere. A variant everyone is diploid at would become one nobody can carry
    once the diploid calls are made nullploid, which is a different test
    """
    _, records = chrx_records
    men: list[str] = []
    kept = []
    for rec in records:
        hemizygous = [s for s, call in rec.samples.items() if len(call["GT"]) == 1]
        if not hemizygous:
            continue
        men = men or hemizygous
        if hemizygous == men:
            kept.append(rec)

    assert len(kept) > 5, "not enough hemizygous chrX variants"
    return men, kept


@pytest.fixture(scope="session")
def chry_bgen_target(
    chrx_records, chry_records, tmp_path_factory
) -> tuple[pathlib.Path, ...]:
    """A mixed-sex chrY file, which only bgen can express: the hemizygous chrX calls
    become chrY calls and the diploid ones become nullploid, because a woman has no
    chrY locus to call. The first man's call at the first variant is dropped so that
    imputation runs at a variant where the women are all-NaN too
    """
    header, _ = chrx_records
    _, records = chry_records
    out_dir = tmp_path_factory.mktemp("chry_bgen")
    path = out_dir / "chrY.bgen"
    samples = list(header.samples)

    with BgenWriter(path, n_samples=len(samples), samples=samples) as bfile:
        for i, rec in enumerate(records):
            probs, ploidy = _probabilities(rec, len(samples))
            probs = probs[:, :2]  # a haploid biallelic call has two probabilities
            probs[ploidy == 2] = np.nan
            ploidy[ploidy == 2] = 0
            if i == 0:
                probs[int(np.flatnonzero(ploidy == 1)[0])] = np.nan
            bfile.add_variant(
                varid=f"var{i}",
                rsid=f"rs{i}",
                chrom="Y",
                pos=rec.pos,
                alleles=[rec.ref, rec.alts[0]],
                genotypes=probs,
                ploidy=ploidy,
                phased=False,
                bit_depth=8,
            )

    subprocess.run(["bgenix", "-index", "-clobber", "-g", str(path)], check=True)
    sample_path = out_dir / "chrY.sample"
    sample_path.write_text(
        "ID_1 ID_2 missing\n0 0 0\n" + "".join(f"{x} {x} 0\n" for x in samples)
    )
    return path, sample_path


def _write_scorefile(path, chrom, records) -> pathlib.Path:
    """Effect allele = ALT and weight = 1, so a score is a sum of dosages"""
    with gzip.open(path, "wt") as f:
        f.write(
            "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\t"
            "effect_type\tis_duplicated\taccession\trow_nr\n"
        )
        for row_nr, rec in enumerate(records):
            f.write(
                f"{chrom}\t{rec.pos}\t{rec.alts[0]}\t{rec.ref}\t1.0\tadditive\tFalse\t"
                f"PGS999993\t{row_nr}\n"
            )
    return path


@pytest.fixture(scope="session")
def chrx_scorefile(chrx_records, tmp_path_factory) -> pathlib.Path:
    _, records = chrx_records
    path = tmp_path_factory.mktemp("chrx_score") / "PGS999993.txt.gz"
    return _write_scorefile(path, "X", records)


@pytest.fixture(scope="session")
def chry_scorefile(chry_records, tmp_path_factory) -> pathlib.Path:
    _, records = chry_records
    path = tmp_path_factory.mktemp("chry_score") / "PGS999993.txt.gz"
    return _write_scorefile(path, "Y", records)


def _score_shard(zip_path, work_dir, scorefile, monkeypatch) -> pl.DataFrame:
    """MATCH -> MATCH_COMBINE -> SCORE -> SCORE_COMBINE over a single shard"""
    limits = ["--max_memory_gb", "16", "--threads", "1"]
    match, combine, score, final = (
        work_dir / name for name in ("match", "match_combine", "score", "combine")
    )
    run_cli(
        [
            "pgsc_calc",
            "match",
            "--zarr_zip_file",
            str(zip_path),
            "--score_paths",
            str(scorefile),
            "--out_dir",
            str(match),
            *limits,
        ],
        monkeypatch,
    )
    run_cli(
        [
            "pgsc_calc",
            "match-combine",
            "--matches",
            str(match / "matches.parquet"),
            "--samples",
            str(match / "samples.parquet"),
            "--score_paths",
            str(scorefile),
            "--sampleset",
            SAMPLESET,
            "--min_overlap",
            "0.75",
            "--out_dir",
            str(combine),
            *limits,
        ],
        monkeypatch,
    )
    run_cli(
        [
            "pgsc_calc",
            "score",
            "--zarr_zip_file",
            str(zip_path),
            "--score_paths",
            str(scorefile),
            "--match_table",
            str(combine / "match_table.parquet"),
            "--out_dir",
            str(score),
            *limits,
        ],
        monkeypatch,
    )
    run_cli(
        [
            "pgsc_calc",
            "combine",
            "--partial_scores",
            str(score / "partial_scores.parquet"),
            "--out_dir",
            str(final),
            *limits,
        ],
        monkeypatch,
    )
    return pl.read_csv(final / "scores.txt.gz").sort("sample_id")


def test_mixed_sex_chrx_bgen_scores_like_the_vcf(
    chrx_vcf_target,
    chrx_bgen_target,
    chrx_scorefile,
    sex_map,
    tmp_path,
    monkeypatch,
):
    """Ploidy arrives on a different channel in each format, and the cache is
    format-blind, so the scores have to be the same number"""
    bgen_path, sample_path = chrx_bgen_target

    run_cli(
        [
            "pgsc_calc",
            "load",
            "--target_genomes",
            f"{bgen_path}:X",
            "--sampleset",
            SAMPLESET,
            "--format",
            "bgen",
            "--cache_dir",
            str(tmp_path / "bgen_cache"),
            "--score_paths",
            str(chrx_scorefile),
            "--bgen_sample_file",
            str(sample_path),
            "--sex_map",
            str(sex_map),
        ],
        monkeypatch,
    )
    bgen_scores = _score_shard(
        tmp_path / "bgen_cache" / "genotypes.zarr.zip",
        tmp_path / "bgen",
        chrx_scorefile,
        monkeypatch,
    )

    vcf_zip = load_target(
        f"{chrx_vcf_target}:X",
        tmp_path / "vcf_cache",
        chrx_scorefile,
        monkeypatch,
        sex_map=sex_map,
    )
    vcf_scores = _score_shard(vcf_zip, tmp_path / "vcf", chrx_scorefile, monkeypatch)

    # men are hemizygous outside the pseudoautosomal regions, so the fixture really
    # does mix ploidy - otherwise this would pass without testing anything
    assert vcf_scores["allele_count"].n_unique() > 1
    assert_frame_equal(bgen_scores, vcf_scores)


def test_a_mixed_sex_chry_bgen_scores(
    chry_records, chry_bgen_target, chry_scorefile, sex_map, tmp_path, monkeypatch
):
    """A woman's chrY call is nullploid rather than missing, which VCF can't say and
    bgen can (ADR 0007): she carries no alleles there, so she takes no dosage, no
    allele count, and no imputed value at a variant where a man's call went missing
    """
    men, records = chry_records
    bgen_path, sample_path = chry_bgen_target

    run_cli(
        [
            "pgsc_calc",
            "load",
            "--target_genomes",
            f"{bgen_path}:Y",
            "--sampleset",
            SAMPLESET,
            "--format",
            "bgen",
            "--cache_dir",
            str(tmp_path / "cache"),
            "--score_paths",
            str(chry_scorefile),
            "--bgen_sample_file",
            str(sample_path),
            "--sex_map",
            str(sex_map),
        ],
        monkeypatch,
    )
    scores = _score_shard(
        tmp_path / "cache" / "genotypes.zarr.zip",
        tmp_path / "score",
        chry_scorefile,
        monkeypatch,
    )

    counts = dict(zip(scores["sample_id"], scores["allele_count"], strict=True))
    dosages = dict(zip(scores["sample_id"], scores["dosage_sum"], strict=True))
    women = set(counts) - set(men)
    assert women and len(men) > 50  # or nothing below is a mixed-sex file

    no_call, called = men[0], men[1:]
    alt = {s: sum(sum(rec.samples[s]["GT"]) for rec in records) for s in men}
    # the effect allele frequency the dropped call is imputed from, taken over the
    # men alone because a woman contributes nothing to either side of it
    freq = sum(sum(records[0].samples[s]["GT"]) for s in called) / len(called)
    assert freq > 0, "an imputed dosage of zero would prove nothing"

    for sample in called:
        assert counts[sample] == len(records)
        assert dosages[sample] == pytest.approx(alt[sample])

    # a missing call is imputed from its own ploidy, and counts no alleles
    assert counts[no_call] == len(records) - 1
    expected = alt[no_call] - sum(records[0].samples[no_call]["GT"]) + freq
    assert dosages[no_call] == pytest.approx(expected)

    # nullploid, so ploidy * freq = 0 at every variant, imputed or not
    assert {counts[s] for s in women} == {0}
    assert {dosages[s] for s in women} == {0}

    # and an average weighted score over no alleles is undefined, not zero
    avg = dict(zip(scores["sample_id"], scores["score_avg"], strict=True))
    assert {avg[s] for s in women} == {None}
