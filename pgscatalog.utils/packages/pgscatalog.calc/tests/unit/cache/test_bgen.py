import pathlib
import shutil

import numpy as np
import pytest

from pgscatalog.calc.lib.constants import MISSING_GENOTYPE_SENTINEL_VALUE
from pgscatalog.calc.lib.cache.genomefiletypes import GenomeFileType
from pgscatalog.calc.lib.cache.targetgenome import TargetGenome
from pgscatalog.calc.lib.cache._bgen import (
    bgen_buffer_variants,
    add_chrom_padding,
    add_chrom_prefix,
    normalise_chrom,
)


#  Chromosome value normalisation
def test_add_chrom_padding_numeric():
    assert add_chrom_padding("1") == "01"
    assert add_chrom_padding(2) == "02"


def test_add_chrom_prefix():
    assert add_chrom_prefix(1) == "chr1"
    assert add_chrom_prefix("X") == "chrX"


def test_normalise_chrom():
    assert normalise_chrom("01") == "1"
    assert normalise_chrom("chr01") == "1"
    assert normalise_chrom("chrX") == "X"


# BGEN buffer variant reading tests


def test_phased_bgen(
    test_positions,
    phased_bgen_path,
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
        idx_path=phased_bgen_path.with_suffix(".bgen.bgi"),
        sample_path=bgen_sample,
        sampleset="test",
    )

    variants.variant_ids
    assert len(variants) == len(test_positions)
    assert len(variants.samples) == n_samples
    assert variants.genotypes.shape == (len(variants), 3202, 2)
    assert np.ptp(variants.genotypes) == np.uint8(1)


def test_unphased_bgen(
    test_positions,
    unphased_bgen_path,
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
        idx_path=unphased_bgen_path.with_suffix(".bgen.bgi"),
        sampleset=bgen_sample,
        sample_path="test",
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
    phased_bgen_path,
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
        idx_path=unphased_bgen_path.with_suffix(".bgen.bgi"),
        sampleset=bgen_samples,
        sample_path="unphased",
    )

    phased = bgen_buffer_variants(
        position_batch=test_positions,
        target_path=phased_bgen_path,
        cache_dir=cache_dir,
        target_chroms=chroms,
        idx_path=phased_bgen_path.with_suffix(".bgen.bgi"),
        sampleset=bgen_samples,
        sample_path="phased",
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
    test_positions, unphased_bgen_path, bgen_sample, tmp_path_factory
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
        idx_path=unphased_bgen_path.with_suffix(".bgen.bgi"),
        sampleset=bgen_sample,
        sample_path="test",
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
    shutil.copy(phased_bgen_path, outdir)
    chroms = [x[0] for x in test_positions]

    new_bgen_path = pathlib.Path(outdir) / phased_bgen_path
    # doesn't exist
    new_bgen_idx = new_bgen_path.with_suffix(".bgen.bgi")

    with pytest.raises(ValueError) as e:
        _ = list(
            bgen_buffer_variants(
                target_path=new_bgen_path,
                position_batch=test_positions,
                cache_dir=cache_dir,
                target_chroms=chroms,
                idx_path=new_bgen_idx.with_suffix(".bgen.bgi"),
                sampleset=bgen_sample,
                sample_path="test",
            )
        )
