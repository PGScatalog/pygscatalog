"""Acceptance tests for imputing missing dosages across mixed ploidy (ADR 0009).

The frequency is pooled across ploidy - total effect alleles over total called
alleles - and the imputed dosage is per (variant, sample): the call's own ploidy
times that one frequency.

Two of the criteria can't be reached through a VCF, because `vcf_get_genotypes`
stores every wholly missing call as diploid whatever its GT length (`.` is both a
haploid and a diploid no-call, and pysam writes `.` for either). So the mixed-sex
chrX file below proves the pooling end to end, and the hemizygous and nullploid
imputed values are driven through `write_dosage_arrays` on a shard built by hand.
"""

import gzip
import pathlib

import numpy as np
import polars as pl
import pysam
import pytest
import zarr
import zarr.storage

from pgscatalog.calc.lib.constants import (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    NO_ALLELE_SENTINEL_VALUE,
)
from pgscatalog.calc.lib.score._dosage import write_dosage_arrays

from .conftest import _is_ambiguous, _write_vcf, run_pipeline

# how many samples lose their call at the imputed variant, which has to clear
# nothing in particular - the sample floor is about the shard, not the variant
N_MISSING = 40


def _write_scorefile(path: pathlib.Path, records: list) -> None:
    """Effect allele = ALT, weight = 1, so score is the sum of effect allele dosages"""
    with gzip.open(path, "wt") as f:
        f.write(
            "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\t"
            "effect_type\tis_duplicated\taccession\trow_nr\n"
        )
        for row_nr, rec in enumerate(records):
            f.write(
                f"X\t{rec.pos}\t{rec.alts[0]}\t{rec.ref}\t1.0\tadditive\tFalse\t"
                f"PGS999993\t{row_nr}\n"
            )


@pytest.fixture(scope="session")
def imputed_chrx(vcf_path, tmp_path_factory):
    """A mixed-sex chrX target where one variant has missing calls.

    Before this ticket LOAD refused the file: a mean dosage over a mixture of
    hemizygous and diploid calls was wrong with nothing to notice it.
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
        ]

    assert len(records) > 10, "not enough chrX test variants"

    # hemizygous and diploid samples both lose their call, so the variant mixes
    # ploidies among the calls that remain
    imputed = records[0]
    calls = list(imputed.samples.items())
    hemizygous = [name for name, call in calls if len(call["GT"]) == 1]
    diploid = [name for name, call in calls if len(call["GT"]) == 2]
    missing = hemizygous[: N_MISSING // 2] + diploid[: N_MISSING // 2]

    # what the file still says at the imputed variant, read before it is blanked
    n_effect_alleles = 0
    n_called_alleles = 0
    for name, call in calls:
        if name in missing:
            continue
        n_effect_alleles += sum(allele == 1 for allele in call["GT"])
        n_called_alleles += len(call["GT"])
    n_called_genotypes = len(calls) - len(missing)

    for name in missing:
        imputed.samples[name]["GT"] = None

    # what each sample carries at every *other* variant, all of them observed
    observed_dosage = {name: 0 for name, _ in calls}
    for rec in records[1:]:
        for name, call in rec.samples.items():
            observed_dosage[name] += sum(allele == 1 for allele in call["GT"])

    directory = tmp_path_factory.mktemp("pooled")
    target = directory / "chrX.vcf.gz"
    _write_vcf(target, header, records)
    scorefile = directory / "PGS999993.txt.gz"
    _write_scorefile(scorefile, records)

    return {
        "target": target,
        "scorefile": scorefile,
        "records": records,
        "missing": missing,
        "observed_dosage": observed_dosage,
        "n_effect_alleles": n_effect_alleles,
        "n_called_alleles": n_called_alleles,
        "n_called_genotypes": n_called_genotypes,
    }


def test_a_mixed_sex_chrx_vcf_with_missing_calls_scores(
    imputed_chrx, tmp_path, monkeypatch
):
    """The frequency isn't halved by the hemizygous samples: it is computed here by
    hand over the same calls, not read back out of the implementation"""
    scores = run_pipeline(
        [(imputed_chrx["target"], "X")],
        tmp_path,
        imputed_chrx["scorefile"],
        monkeypatch,
    )

    # by hand: total effect alleles over total called alleles
    freq = imputed_chrx["n_effect_alleles"] / imputed_chrx["n_called_alleles"]

    # what a ploidy-blind denominator would have produced. Half the callset is
    # hemizygous here, so it is a long way off
    halved = imputed_chrx["n_effect_alleles"] / (imputed_chrx["n_called_genotypes"] * 2)
    assert halved < 0.9 * freq

    imputed_scores = scores.filter(pl.col("sample_id").is_in(imputed_chrx["missing"]))
    assert imputed_scores.shape[0] == N_MISSING

    # a sample's dosage sum is what it was observed to carry everywhere else, plus
    # its imputed call. VCF stores a wholly missing call as diploid whatever its GT
    # length, so every imputed call here is worth two expected copies
    for row in imputed_scores.iter_rows(named=True):
        expected = imputed_chrx["observed_dosage"][row["sample_id"]] + 2 * freq
        assert row["dosage_sum"] == pytest.approx(expected, abs=1e-5), row["sample_id"]
        assert row["dosage_sum"] != pytest.approx(expected - 2 * freq + 2 * halved)

    # a missing call carries no alleles, so it stays out of the denominator of the
    # average weighted score even though it contributes an imputed dosage. Every
    # sample is hemizygous or diploid at every variant in this file
    n = len(imputed_chrx["records"])
    observed = scores.filter(~pl.col("sample_id").is_in(imputed_chrx["missing"]))
    assert set(observed["allele_count"]) == {n, 2 * n}
    assert set(imputed_scores["allele_count"]) == {n - 1, 2 * (n - 1)}


def _score_one_variant(gts):
    """Drive one shard through the real dosage path and read back its dosages"""
    store = zarr.storage.MemoryStore()
    root = zarr.open_group(store=store)
    leaf = root.create_group("test/X/shard.vcf.gz")
    genotypes = leaf.create_array("genotypes", shape=gts.shape, dtype=np.uint8)
    genotypes[:] = gts

    shape = gts.shape[:2]
    pgs = zarr.open_group(store=zarr.storage.MemoryStore())
    dosage = pgs.create_array("dosage", shape=shape, dtype=np.float64)
    missing = pgs.create_array("missing", shape=shape, dtype=np.bool_)
    ploidy = pgs.create_array("ploidy", shape=shape, dtype=np.uint8)

    write_dosage_arrays(
        genotypes_group=root,
        zarr_group="test/X/shard.vcf.gz",
        df=pl.DataFrame(
            {
                "target_row_nr": [0],
                "effect_allele_idx": [1],
                "is_recessive": [False],
                "is_dominant": [False],
            }
        ),
        n_minimum_impute=1,
        dosage_array=dosage,
        missing_array=missing,
        ploidy_array=ploidy,
    )
    return dosage[:], ploidy[:]


def test_a_missing_call_imputes_from_its_own_ploidy():
    """A man's missing chrX call is one expected copy of the effect allele, a
    woman's two - the criterion a VCF can't express, but bgen can"""
    hemi_alt = [1, NO_ALLELE_SENTINEL_VALUE]
    hemi_ref = [0, NO_ALLELE_SENTINEL_VALUE]
    hemi_missing = [MISSING_GENOTYPE_SENTINEL_VALUE, NO_ALLELE_SENTINEL_VALUE]
    diploid_missing = [MISSING_GENOTYPE_SENTINEL_VALUE] * 2
    gts = np.array(
        [[hemi_alt, hemi_ref, [1, 1], [0, 1], hemi_missing, diploid_missing]],
        dtype=np.uint8,
    )

    dosage, ploidy = _score_one_variant(gts)

    # by hand: effect alleles 1 + 0 + 2 + 1 = 4, called alleles 1 + 1 + 2 + 2 = 6
    freq = 4 / 6
    np.testing.assert_array_equal(ploidy[0], [1, 1, 2, 2, 1, 2])
    np.testing.assert_allclose(dosage[0], [1, 0, 2, 1, freq, 2 * freq])


def test_a_nullploid_sample_is_never_given_a_genotype():
    """A woman at a chrY variant contributes nothing to the frequency and receives
    an imputed dosage of zero (ADR 0009)"""
    nullploid = [NO_ALLELE_SENTINEL_VALUE] * 2
    hemi_missing = [MISSING_GENOTYPE_SENTINEL_VALUE, NO_ALLELE_SENTINEL_VALUE]
    gts = np.array(
        [
            [
                [1, NO_ALLELE_SENTINEL_VALUE],
                [0, NO_ALLELE_SENTINEL_VALUE],
                hemi_missing,
                nullploid,
                nullploid,
            ]
        ],
        dtype=np.uint8,
    )

    dosage, ploidy = _score_one_variant(gts)

    # by hand: 1 effect allele over 2 called ones. The nullploid samples are in
    # neither the numerator nor the denominator
    np.testing.assert_array_equal(ploidy[0], [1, 1, 1, 0, 0])
    np.testing.assert_allclose(dosage[0], [1, 0, 0.5, 0, 0])


def test_a_variant_nobody_called_still_stops_the_pipeline():
    """The post-write NaN check catches a frequency taken from no alleles, and says
    which kind of nothing it found"""
    gts = np.full((1, 4, 2), MISSING_GENOTYPE_SENTINEL_VALUE, dtype=np.uint8)

    with pytest.raises(ValueError) as e:
        _score_one_variant(gts)

    assert "no called alleles because every call is missing" in str(e.value)
