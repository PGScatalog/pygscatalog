"""Acceptance tests for MT, the deliberate exception to observed ploidy (ADR 0007).

Everywhere else the caller is the authority on ploidy and the file is honoured.
MT is the one place that gives the wrong answer: callers write the same homoplasmic
variant as `1` or as `1/1`, so honouring the file would let a formatting convention
change a score and two callsets of the same people would disagree. Handlers collapse
the diploid coding instead, and refuse a heterozygous call - heteroplasmy, which has
no defined dosage in a weighted sum.

1000G ships no MT genotypes, so the fixture is real chromosome 1 SNPs rewritten as
homoplasmic MT calls with the contig renamed. What the tests need from it is its
shape: one true dosage per sample, written four different ways.
"""

import gzip
import pathlib
import subprocess

import numpy as np
import pysam
import pytest
from bgen import BgenWriter
from polars.testing import assert_frame_equal

from .conftest import _is_ambiguous, _write_vcf, load_targets, run_serial

# enough MT variants to be worth a sum, few enough to keep four pipeline runs quick
N_MT_VARIANTS = 10


def _write_mt_vcf(path, header, records, dosages, diploid) -> pathlib.Path:
    """Write the same homoplasmic calls as `1` or as `1/1`, on a renamed contig"""
    for rec, dosage in zip(records, dosages, strict=True):
        for name, call in rec.samples.items():
            allele = dosage[name]
            call["GT"] = (allele, allele) if diploid else (allele,)

    raw = path.parent / f"raw_{path.name}"
    _write_vcf(raw, header, records)

    rename = path.parent / f"rename_{path.name}.txt"
    rename.write_text("1\tMT\n")
    pysam.bcftools.annotate(
        "--rename-chrs",
        str(rename),
        "-Oz",
        "-o",
        str(path),
        str(raw),
        catch_stdout=False,
    )
    pysam.bcftools.index(str(path))
    return path


def _write_mt_bgen(out_dir, records, dosages, samples, ploidy) -> tuple:
    """The same calls again, as a bgen recording ploidy 1 or ploidy 2 explicitly"""
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "chrMT.bgen"

    with BgenWriter(path, n_samples=len(samples), samples=samples) as bfile:
        for i, (rec, dosage) in enumerate(zip(records, dosages, strict=True)):
            # unphased: the column is the ALT allele count, so a homoplasmic call is
            # column 1 at ploidy 1 and column 2 at ploidy 2
            probs = np.zeros((len(samples), ploidy + 1))
            for j, name in enumerate(samples):
                probs[j, dosage[name] * ploidy] = 1.0
            bfile.add_variant(
                varid=f"var{i}",
                rsid=f"rs{i}",
                chrom="MT",
                pos=rec.pos,
                alleles=[rec.ref, rec.alts[0]],
                genotypes=probs,
                ploidy=np.full(len(samples), ploidy, dtype=np.uint8),
                phased=False,
                bit_depth=8,
            )

    # BgenWriter leaves its own .bgi behind, but bgenix won't query one it didn't
    # write, so overwrite it
    subprocess.run(["bgenix", "-index", "-clobber", "-g", str(path)], check=True)
    sample_path = out_dir / "chrMT.sample"
    sample_path.write_text(
        "ID_1 ID_2 missing\n0 0 0\n" + "".join(f"{x} {x} 0\n" for x in samples)
    )
    return path, sample_path


@pytest.fixture(scope="session")
def mt_calls(vcf_path):
    """Homoplasmic MT calls: one ALT dosage per sample per variant, which is the
    single fact every coding below has to encode identically"""
    with pysam.VariantFile(str(vcf_path)) as vcf:
        header = vcf.header.copy()
        records = []
        for rec in vcf.fetch("1"):
            if len(rec.ref) != 1 or len(rec.alts) != 1 or len(rec.alts[0]) != 1:
                continue
            if _is_ambiguous(rec.ref, rec.alts[0]):
                continue
            if not all(None not in call["GT"] for call in rec.samples.values()):
                continue
            records.append(rec)
            if len(records) == N_MT_VARIANTS:
                break

    assert len(records) == N_MT_VARIANTS, "not enough MT test variants"

    dosages = [
        {s: int(any(a == 1 for a in call["GT"])) for s, call in rec.samples.items()}
        for rec in records
    ]
    return header, records, dosages


@pytest.fixture(scope="session")
def mt_targets(mt_calls, tmp_path_factory) -> dict:
    """The same calls in four codings, two per format"""
    header, records, dosages = mt_calls
    out_dir = tmp_path_factory.mktemp("mt")
    samples = list(header.samples)

    # the VCF writers mutate GT, so write them in order; the bgen writers read the
    # dosages rather than the records' calls, so they're unaffected either way
    hemi_vcf = _write_mt_vcf(out_dir / "hemi.vcf.gz", header, records, dosages, False)
    diploid_vcf = _write_mt_vcf(out_dir / "dip.vcf.gz", header, records, dosages, True)
    hemi_bgen = _write_mt_bgen(out_dir / "hemi_bgen", records, dosages, samples, 1)
    diploid_bgen = _write_mt_bgen(out_dir / "dip_bgen", records, dosages, samples, 2)

    return {
        "vcf `1`": (hemi_vcf, "vcf", None),
        "vcf `1/1`": (diploid_vcf, "vcf", None),
        "bgen ploidy 1": (hemi_bgen[0], "bgen", hemi_bgen[1]),
        "bgen ploidy 2": (diploid_bgen[0], "bgen", diploid_bgen[1]),
    }


@pytest.fixture(scope="session")
def mt_scorefile(mt_calls, tmp_path_factory) -> pathlib.Path:
    """Effect allele = ALT and weight = 1, so a score is a sum of dosages"""
    _, records, _ = mt_calls
    path = tmp_path_factory.mktemp("mt_score") / "PGS999992.txt.gz"
    with gzip.open(path, "wt") as f:
        f.write(
            "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\t"
            "effect_type\tis_duplicated\taccession\trow_nr\n"
        )
        for row_nr, rec in enumerate(records):
            f.write(
                f"MT\t{rec.pos}\t{rec.alts[0]}\t{rec.ref}\t1.0\tadditive\tFalse\t"
                f"PGS999992\t{row_nr}\n"
            )
    return path


def test_the_same_homoplasmic_mt_variant_scores_the_same_however_its_written(
    mt_calls, mt_targets, mt_scorefile, tmp_path, monkeypatch
):
    """The point of ticket 0023: `1` and `1/1` are one caller's formatting choice
    and the other's, and both have to produce the same dosage - in both formats,
    where ploidy arrives on entirely different channels.

    MT allele counts come out at 1 per matched variant, not 2, because the collapsed
    call is hemizygous and ploidy is derived from the stored array.
    """
    _, _, dosages = mt_calls
    expected = {
        sample: sum(dosage[sample] for dosage in dosages) for sample in dosages[0]
    }

    scores = {}
    for i, (label, (path, fmt, sample_file)) in enumerate(mt_targets.items()):
        scores[label] = run_serial(
            [(path, "MT")],
            tmp_path / f"run{i}",
            mt_scorefile,
            monkeypatch,
            fmt=fmt,
            sample_file=sample_file,
        ).sort("sample_id")

    # the fixture would pass vacuously if every sample carried the same dosage
    assert len(set(expected.values())) > 1

    reference = scores["vcf `1`"]
    for label, frame in scores.items():
        assert_frame_equal(frame, reference), label

    assert reference.shape[0] == 3202
    assert reference["n_matched"].unique().to_list() == [N_MT_VARIANTS]
    # one allele per matched MT variant, not two: mitochondria are hemizygous
    assert reference["allele_count"].unique().to_list() == [N_MT_VARIANTS]

    for row in reference.iter_rows(named=True):
        assert row["dosage_sum"] == expected[row["sample_id"]], row["sample_id"]


def test_a_heterozygous_mt_call_is_refused_by_name(
    mt_targets, mt_scorefile, tmp_path, monkeypatch
):
    """A het MT call is heteroplasmy - a within-sample mixture - and polygenic
    scoring has no defined dosage for it, so LOAD refuses rather than picking a
    number or calling it missing and imputing one (ADR 0007)"""
    path, _, _ = mt_targets["vcf `1/1`"]
    with pysam.VariantFile(str(path)) as vcf:
        header = vcf.header.copy()
        records = list(vcf.fetch("MT"))

    list(records[0].samples.values())[0]["GT"] = (0, 1)
    het_path = tmp_path / "heteroplasmic.vcf.gz"
    _write_vcf(het_path, header, records)

    with pytest.raises(ValueError) as excinfo:
        load_targets([(het_path, "MT")], tmp_path / "cache", mt_scorefile, monkeypatch)

    message = str(excinfo.value)
    # a user reading this has to see the pipeline refusing, not failing
    assert "heteroplasmy" in message
    assert "MT:" in message
