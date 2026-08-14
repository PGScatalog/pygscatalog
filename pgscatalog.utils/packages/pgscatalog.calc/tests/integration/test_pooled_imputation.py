"""Acceptance tests for imputing missing dosages across mixed ploidy (ADR 0009).

The frequency is pooled across ploidy - total effect alleles over total called
alleles - and the imputed dosage is per (variant, sample): the call's own ploidy
times that one frequency.

The reader now reads a wholly missing call's ploidy from its GT length (ADR 0007),
so a hemizygous no-call reaches the pipeline through a real VCF: the mixed-sex chrX
file below blanks men and women alike at one variant, and each is imputed at the
ploidy its own `GT` declared.

Nullploidy still has no VCF spelling at all, so that one imputed value is driven
through `write_dosage_arrays` on a shard built by hand.
"""

import gzip
import pathlib

import numpy as np
import polars as pl
import pysam
import pytest
import zarr
import zarr.storage

from pgscatalog.calc.cli.load import unzip_zarr
from pgscatalog.calc.lib.constants import (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    NO_ALLELE_SENTINEL_VALUE,
)
from pgscatalog.calc.lib.score._dosage import write_dosage_arrays

from .conftest import (
    SAMPLESET,
    _is_ambiguous,
    _subset_samples,
    _write_vcf,
    load_target,
    load_targets,
    run_serial,
)

# how many samples lose their call at the imputed variant, which has to clear
# nothing in particular - the sample floor is about the shard, not the variant
N_MISSING = 40

# enough chrX variants to make a readable dosage sum, few enough to stay quick
N_HEMIZYGOUS_VARIANTS = 10


def _write_scorefile(path: pathlib.Path, records: list, chrom: str = "X") -> None:
    """Effect allele = ALT, weight = 1, so score is the sum of effect allele dosages"""
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


def _diploid_records(vcf_path, chrom: str, n: int) -> tuple:
    """Autosomal records, where every call in the file is diploid and a `.` has
    nothing at the variant to explain it"""
    with pysam.VariantFile(str(vcf_path)) as vcf:
        header = vcf.header.copy()
        records = []
        for rec in vcf.fetch(chrom):
            if len(rec.ref) != 1 or len(rec.alts) != 1 or len(rec.alts[0]) != 1:
                continue
            if _is_ambiguous(rec.ref, rec.alts[0]):
                continue
            records.append(rec)
            if len(records) == n:
                break

    assert len(records) == n, "not enough diploid test variants"
    assert all(len(call["GT"]) == 2 for call in records[0].samples.values())
    return header, records


def _hemizygous_chrx_records(vcf_path, n: int) -> tuple:
    """Non-pseudoautosomal chrX records where the same samples are hemizygous in
    every one of them, plus that sample list"""
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
            if len(records) == n:
                break

    assert len(records) == n, "not enough hemizygous test variants"
    return header, records, sorted(males)


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
    # ploidies among the calls that remain. Each blanked call is written the way a
    # conformant writer would write it: `.` for a man, `./.` for a woman
    imputed = records[0]
    calls = list(imputed.samples.items())
    hemizygous = [name for name, call in calls if len(call["GT"]) == 1]
    diploid = [name for name, call in calls if len(call["GT"]) == 2]
    assert hemizygous and diploid, "not a mixed-sex variant"
    missing = hemizygous[: N_MISSING // 2] + diploid[: N_MISSING // 2]
    ploidy = {name: len(call["GT"]) for name, call in calls}

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
        imputed.samples[name]["GT"] = (None,) * ploidy[name]

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
        "ploidy": ploidy,
        "observed_dosage": observed_dosage,
        "n_effect_alleles": n_effect_alleles,
        "n_called_alleles": n_called_alleles,
        "n_called_genotypes": n_called_genotypes,
    }


def test_a_mixed_sex_chrx_vcf_with_missing_calls_scores(
    imputed_chrx, sex_map, tmp_path, monkeypatch
):
    """The frequency isn't halved by the hemizygous samples: it is computed here by
    hand over the same calls, not read back out of the implementation"""
    scores = run_serial(
        [(imputed_chrx["target"], "X")],
        tmp_path,
        imputed_chrx["scorefile"],
        monkeypatch,
        sex_map=sex_map,
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
    # its imputed call - which is worth its own ploidy in expected copies, one for
    # the man whose call was blanked to `.` and two for the woman blanked to `./.`
    for row in imputed_scores.iter_rows(named=True):
        sample = row["sample_id"]
        ploidy = imputed_chrx["ploidy"][sample]
        expected = imputed_chrx["observed_dosage"][sample] + ploidy * freq
        assert row["dosage_sum"] == pytest.approx(expected, abs=1e-5), sample
        assert row["dosage_sum"] != pytest.approx(
            expected - ploidy * freq + ploidy * halved
        )

    # a missing call carries no alleles, so it stays out of the denominator of the
    # average weighted score even though it contributes an imputed dosage. Every
    # sample is hemizygous or diploid at every variant in this file
    n = len(imputed_chrx["records"])
    observed = scores.filter(~pl.col("sample_id").is_in(imputed_chrx["missing"]))
    assert set(observed["allele_count"]) == {n, 2 * n}
    assert set(imputed_scores["allele_count"]) == {n - 1, 2 * (n - 1)}


@pytest.fixture(scope="session")
def hemizygous_chrx(vcf_path, tmp_path_factory):
    """A male-only chrX target where one variant has missing calls.

    Every call in the file is hemizygous, so a `.` in it is unambiguously a
    hemizygous no-call and is stored as one - the case ADR 0007's encoding table
    calls [255, 254] and the reader used to flatten to a diploid missing call.
    """
    out_dir = tmp_path_factory.mktemp("hemizygous")
    header, records, males = _hemizygous_chrx_records(vcf_path, N_HEMIZYGOUS_VARIANTS)

    everyone = out_dir / "chrX_all.vcf.gz"
    _write_vcf(everyone, header, records)
    male_only = _subset_samples(everyone, out_dir / "male_only.vcf.gz", males)

    with pysam.VariantFile(str(male_only)) as vcf:
        male_header = vcf.header.copy()
        male_records = list(vcf.fetch("X"))

    imputed = male_records[0]
    calls = list(imputed.samples.items())
    assert all(len(call["GT"]) == 1 for _, call in calls), "not a male-only file"
    missing = [name for name, _ in calls][:N_MISSING]

    # what the file still says at the imputed variant, read before it is blanked.
    # Every remaining call is hemizygous, so one call is one allele
    n_effect_alleles = sum(
        sum(allele == 1 for allele in call["GT"])
        for name, call in calls
        if name not in missing
    )
    n_called_alleles = len(calls) - len(missing)

    for name in missing:
        imputed.samples[name]["GT"] = (None,)

    observed_dosage = {name: 0 for name, _ in calls}
    for rec in male_records[1:]:
        for name, call in rec.samples.items():
            observed_dosage[name] += sum(allele == 1 for allele in call["GT"])

    target = out_dir / "chrX.vcf.gz"
    _write_vcf(target, male_header, male_records)
    scorefile = out_dir / "PGS999993.txt.gz"
    _write_scorefile(scorefile, male_records)

    return {
        "target": target,
        "scorefile": scorefile,
        "records": male_records,
        "missing": missing,
        "observed_dosage": observed_dosage,
        "freq": n_effect_alleles / n_called_alleles,
    }


def test_a_mans_missing_chrx_call_imputes_at_ploidy_one(
    hemizygous_chrx, sex_map, tmp_path, monkeypatch
):
    """A hemizygous no-call is one expected copy of the effect allele, not two, and
    it counts one allele toward the sample's denominator - through the real reader,
    not a shard built by hand.

    The male-only file is the unambiguous case: nothing in it is diploid, so there
    is nothing a `.` could be but a hemizygous no-call.
    """
    scores = run_serial(
        [(hemizygous_chrx["target"], "X")],
        tmp_path,
        hemizygous_chrx["scorefile"],
        monkeypatch,
        sex_map=sex_map,
    )

    freq = hemizygous_chrx["freq"]
    missing = hemizygous_chrx["missing"]
    imputed_scores = scores.filter(pl.col("sample_id").is_in(missing))
    assert imputed_scores.shape[0] == N_MISSING

    for row in imputed_scores.iter_rows(named=True):
        expected = hemizygous_chrx["observed_dosage"][row["sample_id"]] + freq
        assert row["dosage_sum"] == pytest.approx(expected, abs=1e-5), row["sample_id"]
        # what storing the `.` as a diploid missing call produced
        assert row["dosage_sum"] != pytest.approx(expected + freq, abs=1e-5)

    # the missing call carries no alleles at all, and every other call carries one
    n = len(hemizygous_chrx["records"])
    observed = scores.filter(~pl.col("sample_id").is_in(hemizygous_chrx["missing"]))
    assert set(observed["allele_count"]) == {n}
    assert set(imputed_scores["allele_count"]) == {n - 1}


def test_a_dot_where_nothing_is_hemizygous_stops_the_load(
    vcf_path, tmp_path, monkeypatch
):
    """The boundary of the trust: on an autosome every call is diploid, so a `.` has
    no reading at all and can only be a writer that dropped an allele. LOAD says so
    instead of imputing those samples at half the dosage they can carry"""
    header, records = _diploid_records(vcf_path, "1", 2)

    blanked = list(records[0].samples)[:5]
    for name in blanked:
        records[0].samples[name]["GT"] = (None,)

    target = tmp_path / "autosome.vcf.gz"
    _write_vcf(target, header, records)
    scorefile = tmp_path / "PGS999993.txt.gz"
    _write_scorefile(scorefile, records, chrom="1")

    with pytest.raises(ValueError) as excinfo:
        load_targets([(target, "1")], tmp_path / "cache", scorefile, monkeypatch)

    message = str(excinfo.value)
    assert f"1:{records[0].pos}" in message
    assert "Non-conformant VCF" in message


def test_a_dot_at_a_mixed_sex_chrx_variant_is_stored_hemizygous(
    vcf_path, sex_map, tmp_path, monkeypatch
):
    """The mirror of it: the same `.` at a variant where some sample really is
    hemizygous is ordinary conformant output, so it loads and is stored [255, 254].

    Refusing this would refuse most real chrX data, since a mixed-sex file with any
    missing call would be unreadable (ADR 0007).
    """
    header, records, males = _hemizygous_chrx_records(vcf_path, 2)

    blanked = males[:5]
    for name in blanked:
        records[0].samples[name]["GT"] = (None,)

    target = tmp_path / "mixed_sex_chrX.vcf.gz"
    _write_vcf(target, header, records)
    scorefile = tmp_path / "PGS999993.txt.gz"
    _write_scorefile(scorefile, records)

    zip_path = load_target(
        f"{target}:X", tmp_path / "cache", scorefile, monkeypatch, sex_map=sex_map
    )
    unzip_zarr(zip_path, tmp_path / "unzipped")
    store = zarr.storage.LocalStore(tmp_path / "unzipped" / "genotypes.zarr")
    root = zarr.open_group(store=store, mode="r")
    leaf = root[f"{SAMPLESET}/X/{target.name}"]

    samples = list(leaf.attrs["samples"])
    positions = leaf["meta/chr_pos"][:]
    row = int(np.flatnonzero(positions == records[0].pos)[0])
    genotypes = leaf["genotypes"][row]

    for name in blanked:
        call = genotypes[samples.index(name)]
        np.testing.assert_array_equal(
            call, [MISSING_GENOTYPE_SENTINEL_VALUE, NO_ALLELE_SENTINEL_VALUE]
        )

    # the women at the same variant kept both their allele slots
    woman = next(s for s in samples if s not in males)
    assert genotypes[samples.index(woman)][1] != NO_ALLELE_SENTINEL_VALUE


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
