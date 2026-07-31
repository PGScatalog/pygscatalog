import shutil
import sqlite3
import struct
import zlib

import numpy as np
import pytest
from bgen import BgenReader, BgenWriter

from pgscatalog.calc.lib.cache._bgen import (
    add_chrom_padding,
    add_chrom_prefix,
    bgen_buffer_variants,
    bgen_to_hardcall,
    parse_target_variants,
)
from pgscatalog.calc.lib.constants import (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    NO_ALLELE_SENTINEL_VALUE,
)

# probability rows for a biallelic site. A haploid call has two probabilities under
# both the phased and the unphased representation, which is what makes haploid
# hard-calling identical either way (ADR 0007)
HAP_REF, HAP_ALT = [1.0, 0.0], [0.0, 1.0]
HOM_REF, HET, HOM_ALT = [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]
PHASED_HET, PHASED_HOM_ALT = [1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, 1.0]
NAN = float("nan")


#  Chromosome value normalisation
def test_add_chrom_padding_numeric():
    assert add_chrom_padding("1") == "01"
    assert add_chrom_padding(2) == "02"


def test_add_chrom_prefix():
    assert add_chrom_prefix(1) == "chr1"
    assert add_chrom_prefix("X") == "chrX"


# BGEN buffer variant reading tests


def test_padded_bgen_loads_unchanged(
    test_positions,
    phased_bgen_path,
    phased_bgen_index_path,
    bgen_sample,
    tmp_path_factory,
):
    """Zero-padded targets (UK Biobank style) load exactly like unpadded ones.

    The query still goes out in the target's dialect, and what comes back is
    canonicalised, so the cache can't tell the two apart (ADR 0008).
    """
    padded_index = shutil.copy(
        phased_bgen_index_path, tmp_path_factory.mktemp("padded") / "padded.bgen.bgi"
    )
    with sqlite3.connect(padded_index) as conn:
        conn.execute("UPDATE Variant SET chromosome = printf('%02d', chromosome)")

    kwargs = {
        "position_batch": test_positions,
        "target_path": phased_bgen_path,
        "sample_path": bgen_sample,
        "sampleset": "test",
    }
    padded = bgen_buffer_variants(
        cache_dir=tmp_path_factory.mktemp("padded_cache"),
        target_chroms=["01", "08", "10"],
        idx_path=padded_index,
        **kwargs,
    )
    unpadded = bgen_buffer_variants(
        cache_dir=tmp_path_factory.mktemp("cache"),
        target_chroms=[x[0] for x in test_positions],
        idx_path=phased_bgen_index_path,
        **kwargs,
    )

    # keyed by variant, because a padded index sorts 01/08/10 where an unpadded one
    # sorts 1/10/8 - the rows are the same, the order they come back in isn't
    padded_gts = dict(zip(padded.variant_ids, padded.genotypes, strict=True))
    unpadded_gts = dict(zip(unpadded.variant_ids, unpadded.genotypes, strict=True))

    assert padded_gts.keys() == unpadded_gts.keys()
    assert all(np.array_equal(padded_gts[k], unpadded_gts[k]) for k in padded_gts)


def test_phased_bgen(
    test_positions,
    phased_bgen_path,
    phased_bgen_index_path,
    bgen_sample,
    n_samples,
    tmp_path_factory,
    sample_ids,
):
    cache_dir = tmp_path_factory.mktemp("cache")
    chroms = [x[0] for x in test_positions]
    variants = bgen_buffer_variants(
        position_batch=test_positions,
        target_path=phased_bgen_path,
        cache_dir=cache_dir,
        target_chroms=chroms,
        idx_path=phased_bgen_index_path,
        sample_path=bgen_sample,
        sampleset="test",
    )

    assert len(variants) == len(test_positions)
    assert len(variants.samples) == n_samples
    assert variants.genotypes.shape == (len(variants), 3202, 2)
    assert np.ptp(variants.genotypes) == np.uint8(1)


def test_unphased_bgen(
    test_positions,
    unphased_bgen_path,
    unphased_bgen_index_path,
    bgen_sample,
    n_samples,
    tmp_path_factory,
    sample_ids,
):
    cache_dir = tmp_path_factory.mktemp("cache")
    chroms = [x[0] for x in test_positions]
    variants = bgen_buffer_variants(
        position_batch=test_positions,
        target_path=unphased_bgen_path,
        cache_dir=cache_dir,
        target_chroms=chroms,
        idx_path=unphased_bgen_index_path,
        sampleset="test",
        sample_path=bgen_sample,
    )

    assert len(variants) == len(test_positions)
    assert len(variants.samples) == n_samples
    assert list(variants.samples) == list(sample_ids)

    gts = variants.genotypes
    assert gts.shape == (3, 3202, 2)
    assert np.ptp(gts) == np.uint8(1)


def test_phase_equality(
    test_positions,
    unphased_bgen_path,
    unphased_bgen_index_path,
    phased_bgen_path,
    phased_bgen_index_path,
    bgen_sample,
    bgen_samples,
    n_samples,
    tmp_path_factory,
):
    """Test that genotypes are called the same for phased and unphased bgens."""
    cache_dir = tmp_path_factory.mktemp("cache")
    chroms = [x[0] for x in test_positions]

    unphased = bgen_buffer_variants(
        position_batch=test_positions,
        target_path=unphased_bgen_path,
        cache_dir=cache_dir,
        target_chroms=chroms,
        idx_path=unphased_bgen_index_path,
        sampleset="unphased",
        sample_path=bgen_sample,
    )

    phased = bgen_buffer_variants(
        position_batch=test_positions,
        target_path=phased_bgen_path,
        cache_dir=cache_dir,
        target_chroms=chroms,
        idx_path=phased_bgen_index_path,
        sampleset="phased",
        sample_path=bgen_sample,
    )

    assert unphased.genotypes.shape == phased.genotypes.shape

    unphased_gts = unphased.genotypes
    phased_gts = phased.genotypes

    assert unphased_gts.shape == (len(test_positions), n_samples, 2)
    assert phased_gts.shape == (len(test_positions), n_samples, 2)

    for unphased, phased in zip(unphased_gts, phased_gts, strict=True):
        # but precise gt order for heterozygous (0, 1) or (1, 0) can differ
        # (because it's inferred in unphased samples)
        assert np.array_equal(np.sum(unphased, axis=1), np.sum(phased, axis=1))


def test_bgen_with_missing(
    test_positions,
    unphased_bgen_path,
    unphased_bgen_index_path,
    bgen_sample,
    tmp_path_factory,
):
    cache_dir = tmp_path_factory.mktemp("cache")
    missing_chrom, missing_pos = "1", 9067156  # definitely missing
    positions = test_positions + [(missing_chrom, missing_pos)]
    chroms = [x[0] for x in positions]

    variants = bgen_buffer_variants(
        position_batch=positions,
        target_path=unphased_bgen_path,
        cache_dir=cache_dir,
        target_chroms=chroms,
        idx_path=unphased_bgen_index_path,
        sampleset="test",
        sample_path=bgen_sample,
    )

    # test that missing variants are yielded properly
    assert len(variants) == len(positions)
    missing_index = positions.index((missing_chrom, missing_pos))
    assert variants.variant_metadata.chr_name[missing_index] == "1"
    assert variants.variant_metadata.ref[missing_index] is None
    assert variants.variant_metadata.alts[missing_index] is None
    # Find the record for the missing variant
    meta = variants.variant_metadata
    chr_arr = np.asarray(meta.chr_name)
    pos_arr = np.asarray(meta.chr_pos)

    mask = (chr_arr == missing_chrom) & (pos_arr == missing_pos)
    idx = int(np.where(mask)[0][0])
    print("Index of missing variant:", idx)
    assert meta.ref[idx] is None
    assert meta.alts[idx] is None
    assert np.all(variants.genotypes[idx] == MISSING_GENOTYPE_SENTINEL_VALUE)


def test_bgen_no_index(test_positions, phased_bgen_path, bgen_sample, tmp_path_factory):
    cache_dir = tmp_path_factory.mktemp("cache")
    outdir = tmp_path_factory.mktemp("bgen_no_index")
    new_bgen_path = shutil.copy(phased_bgen_path, outdir)
    chroms = [x[0] for x in test_positions]

    # doesn't exist
    new_bgen_idx = new_bgen_path + ".bgi"

    with pytest.raises(ValueError):
        _ = bgen_buffer_variants(
            target_path=new_bgen_path,
            position_batch=test_positions,
            cache_dir=cache_dir,
            target_chroms=chroms,
            idx_path=new_bgen_idx,
            sampleset="test",
            sample_path=bgen_sample,
        )


# Ploidy (ADR 0007). These read a bgen directly instead of going through bgenix,
# because the interesting thing is the hard call, not the index query


def _write_sample_file(path, samples):
    """A .sample file with the ID in both identifier columns"""
    path.write_text(
        "ID_1 ID_2 missing\n0 0 0\n" + "".join(f"{x} {x} 0\n" for x in samples)
    )
    return path


def _write_bgen(tmp_path, variants, n_samples=4):
    """Write a bgen holding exactly the calls a test needs.

    ``variants`` are (chrom, pos, probabilities, ploidy, phased) tuples. Ploidy is
    an explicit field in layout 2, so the writer takes it rather than deriving it.
    """
    samples = [f"sample{i}" for i in range(n_samples)]
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "test.bgen"
    with BgenWriter(path, n_samples=n_samples, samples=samples) as bfile:
        for i, (chrom, pos, probs, ploidy, phased) in enumerate(variants):
            bfile.add_variant(
                varid=f"var{i}",
                rsid=f"rs{i}",
                chrom=chrom,
                pos=pos,
                alleles=["A", "G"],
                genotypes=np.array(probs, dtype=np.float64),
                ploidy=np.array(ploidy, dtype=np.uint8),
                phased=phased,
                bit_depth=8,
            )
    return path, _write_sample_file(tmp_path / "test.sample", samples)


def _write_layout1_bgen(path, samples, variants):
    """Hand-roll a layout 1 bgen, which no writer in the ecosystem produces.

    Layout 1 has no per-sample ploidy field at all - three probabilities per sample,
    always - which is the whole reason a non-autosome can't be read from one.
    """
    n = len(samples)
    header = (
        struct.pack("<III", 20, len(variants), n)
        + b"bgen"
        + struct.pack("<I", 1 | (1 << 2))  # zlib compressed, layout 1
    )
    blocks = b""
    for i, (chrom, pos, probs) in enumerate(variants):
        block = struct.pack("<I", n)
        for text in (f"var{i}", f"rs{i}", chrom):
            block += struct.pack("<H", len(text)) + text.encode()
        block += struct.pack("<I", pos)
        for allele in ("A", "G"):
            block += struct.pack("<I", 1) + allele.encode()
        raw = np.round(np.array(probs) * 32768).astype("<u2").tobytes()
        compressed = zlib.compress(raw)
        blocks += block + struct.pack("<I", len(compressed)) + compressed

    path.write_bytes(struct.pack("<I", 20) + header + blocks)
    return path, _write_sample_file(path.with_suffix(".sample"), samples)


def _parse(path, sample_path, positions):
    return parse_target_variants(
        target_path=path,
        buffered_path=path,
        bgen_sample_path=sample_path,
        scoring_file_regions=positions,
        sampleset="test",
    )


def test_hemizygous_call_fills_the_slot_it_doesnt_have(tmp_path):
    """Ploidy comes from BgenVar.ploidy, and a haploid call keeps one allele"""
    path, sample_path = _write_bgen(
        tmp_path,
        [
            (
                "X",
                100,
                [HAP_REF + [NAN], HAP_ALT + [NAN], HET, HOM_ALT],
                [1, 1, 2, 2],
                False,
            )
        ],
    )
    variants = _parse(path, sample_path, [("X", 100)])

    assert np.array_equal(
        variants.genotypes[0],
        [[0, NO_ALLELE_SENTINEL_VALUE], [1, NO_ALLELE_SENTINEL_VALUE], [0, 1], [1, 1]],
    )


def test_phased_and_unphased_haploid_calls_agree(tmp_path):
    """Two probabilities mean the same thing under both representations"""
    calls = [HAP_REF, HAP_ALT, HAP_ALT, HAP_REF]
    unphased, unphased_samples = _write_bgen(
        tmp_path / "unphased", [("X", 100, calls, [1, 1, 1, 1], False)]
    )
    phased, phased_samples = _write_bgen(
        tmp_path / "phased", [("X", 100, calls, [1, 1, 1, 1], True)]
    )

    assert np.array_equal(
        _parse(unphased, unphased_samples, [("X", 100)]).genotypes,
        _parse(phased, phased_samples, [("X", 100)]).genotypes,
    )


def test_phasing_is_read_per_variant(tmp_path):
    """is_phased is a property of a variant, not of a file, so one file can mix
    both and each variant has to pick its own hard-call path"""
    path, sample_path = _write_bgen(
        tmp_path,
        [
            (
                "X",
                100,
                [
                    HAP_REF + [NAN, NAN],
                    HAP_ALT + [NAN, NAN],
                    PHASED_HET,
                    PHASED_HOM_ALT,
                ],
                [1, 1, 2, 2],
                1,
            ),
            (
                "X",
                200,
                [HAP_REF + [NAN], HAP_ALT + [NAN], HET, HOM_ALT],
                [1, 1, 2, 2],
                False,
            ),
        ],
    )
    variants = _parse(path, sample_path, [("X", 100), ("X", 200)])

    expected = [
        [0, NO_ALLELE_SENTINEL_VALUE],
        [1, NO_ALLELE_SENTINEL_VALUE],
        [0, 1],
        [1, 1],
    ]
    assert np.array_equal(variants.genotypes[0], expected)
    assert np.array_equal(variants.genotypes[1], expected)


def test_ploidy_above_two_raises_naming_the_ploidy(tmp_path):
    """The genotype axis is two wide, so a triploid call can't be held. Keeping the
    first two alleles would turn it into a plausible score"""
    triploid = [[0.7, 0.1, 0.1, 0.1]] * 4
    path, sample_path = _write_bgen(
        tmp_path, [("X", 100, triploid, [3, 3, 3, 3], False)]
    )

    with pytest.raises(ValueError, match="Ploidy 3"):
        _parse(path, sample_path, [("X", 100)])


def test_a_missing_call_is_missing_not_nullploid(tmp_path):
    """There's no missingness flag, so a missing call is a wholly NaN row - and a
    missing haploid call keeps its ploidy, because bgen still records it"""
    path, sample_path = _write_bgen(
        tmp_path,
        [
            (
                "X",
                100,
                [[NAN, NAN, NAN], HAP_ALT + [NAN], [NAN, NAN, NAN], HOM_ALT],
                [1, 1, 2, 2],
                0,
            )
        ],
    )
    genotypes = _parse(path, sample_path, [("X", 100)]).genotypes[0]

    assert np.array_equal(
        genotypes,
        [
            [MISSING_GENOTYPE_SENTINEL_VALUE, NO_ALLELE_SENTINEL_VALUE],
            [1, NO_ALLELE_SENTINEL_VALUE],
            [MISSING_GENOTYPE_SENTINEL_VALUE] * 2,
            [1, 1],
        ],
    )
    # ploidy is derived from the array, so a missing call still has one (ADR 0007)
    ploidy = np.sum(genotypes != NO_ALLELE_SENTINEL_VALUE, axis=1)
    assert np.array_equal(ploidy, [1, 1, 2, 2])


def _mixed_sex_chry(tmp_path):
    """A chrY variant with a woman, a man whose call went missing, and two called
    men. The first two rows are both wholly NaN, which is the point"""
    return _write_bgen(
        tmp_path,
        [("Y", 100, [[NAN, NAN], [NAN, NAN], HAP_ALT, HAP_REF], [0, 1, 1, 1], False)],
    )


def test_a_nullploid_call_is_not_a_missing_one(tmp_path):
    """A woman has no chrY locus to call, so she reaches the cache nullploid rather
    than missing - which is what keeps her out of imputation later"""
    path, sample_path = _mixed_sex_chry(tmp_path)
    genotypes = _parse(path, sample_path, [("Y", 100)]).genotypes[0]

    assert np.array_equal(
        genotypes,
        [
            [NO_ALLELE_SENTINEL_VALUE] * 2,
            [MISSING_GENOTYPE_SENTINEL_VALUE, NO_ALLELE_SENTINEL_VALUE],
            [1, NO_ALLELE_SENTINEL_VALUE],
            [0, NO_ALLELE_SENTINEL_VALUE],
        ],
    )
    ploidy = np.sum(genotypes != NO_ALLELE_SENTINEL_VALUE, axis=1)
    assert np.array_equal(ploidy, [0, 1, 1, 1])


def test_nullploid_and_missing_are_told_apart_by_the_ploidy_accessor(tmp_path):
    """The probabilities can't tell them apart and must never be asked to: both
    rows are identically NaN, so BgenVar.ploidy is the only source (ADR 0007)"""
    path, sample_path = _mixed_sex_chry(tmp_path)
    with BgenReader(str(path), sample_path=str(sample_path)) as bfile:
        variant = next(iter(bfile))
        probs, ploidy = variant.probabilities, variant.ploidy

    assert np.isnan(probs[:2]).all()
    np.testing.assert_array_equal(probs[0], probs[1])

    gts = bgen_to_hardcall(probs=probs, ploidy=ploidy, is_phased=False)
    assert not np.array_equal(gts[0], gts[1])

    # swap what the accessor says and the two calls swap with it, with the
    # probabilities untouched
    swapped = bgen_to_hardcall(
        probs=probs, ploidy=np.array([1, 0, 1, 1], np.uint8), is_phased=False
    )
    assert np.array_equal(swapped[0], gts[1])
    assert np.array_equal(swapped[1], gts[0])


def test_diploid_coded_mt_calls_collapse_to_hemizygous(tmp_path):
    """MT is the one place ploidy isn't observed from the file (ADR 0007): the same
    homoplasmic variant is ploidy 1 from one caller and a homozygous ploidy 2 call
    from another, and both have to reach the cache as the same hemizygous call"""
    hemizygous, sample_path = _write_bgen(
        tmp_path / "hemizygous",
        [("MT", 100, [HAP_REF, HAP_ALT, HAP_ALT, HAP_REF], [1, 1, 1, 1], False)],
    )
    diploid, diploid_samples = _write_bgen(
        tmp_path / "diploid",
        [("MT", 100, [HOM_REF, HOM_ALT, HOM_ALT, HOM_REF], [2, 2, 2, 2], False)],
    )

    expected = [
        [0, NO_ALLELE_SENTINEL_VALUE],
        [1, NO_ALLELE_SENTINEL_VALUE],
        [1, NO_ALLELE_SENTINEL_VALUE],
        [0, NO_ALLELE_SENTINEL_VALUE],
    ]
    assert np.array_equal(
        _parse(hemizygous, sample_path, [("MT", 100)]).genotypes[0], expected
    )
    assert np.array_equal(
        _parse(diploid, diploid_samples, [("MT", 100)]).genotypes[0], expected
    )


def test_a_heterozygous_mt_call_names_heteroplasmy(tmp_path):
    """A within-sample mixture has no defined dosage in a weighted sum, so the
    pipeline refuses it rather than picking a number (ADR 0007)"""
    path, sample_path = _write_bgen(
        tmp_path, [("MT", 100, [HOM_REF, HET, HOM_ALT, HOM_ALT], [2, 2, 2, 2], False)]
    )

    with pytest.raises(ValueError, match="heteroplasmy"):
        _parse(path, sample_path, [("MT", 100)])


def test_layout_1_bgen_rejects_a_non_autosome(tmp_path):
    path, sample_path = _write_layout1_bgen(
        tmp_path / "layout1.bgen",
        ["sample0", "sample1"],
        [("X", 100, [HOM_REF, HET])],
    )

    with pytest.raises(ValueError, match="layout 1"):
        _parse(path, sample_path, [("X", 100)])


def test_layout_1_bgen_still_reads_an_autosome(tmp_path):
    """Diploid by assumption is the right answer on an autosome"""
    path, sample_path = _write_layout1_bgen(
        tmp_path / "layout1.bgen",
        ["sample0", "sample1"],
        [("1", 100, [HOM_REF, HET])],
    )

    variants = _parse(path, sample_path, [("1", 100)])

    assert np.array_equal(variants.genotypes[0], [[0, 0], [0, 1]])
