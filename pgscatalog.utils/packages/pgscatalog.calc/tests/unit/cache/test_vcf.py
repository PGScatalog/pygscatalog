import shutil
import pathlib

import pysam
import pytest
import numpy as np

from pgscatalog.calc.lib.cache.targetvariants import TargetVariants, add_missing_positions_to_lists
from pgscatalog.calc.lib.constants import (
    MISSING_GENOTYPE_TUPLE,
)

from pgscatalog.calc.lib.cache._vcf import (
    vcf_get_genotypes, 
    vcf_get_sample_list,
    vcf_buffer_variants,
    parse_target_variants,
    run_bcftools_view
)

@pytest.fixture
def test_positions() -> list[tuple[str, int]]:
    """Some positions in tiny1000G"""
    return [("8", 116938529), ("1", 44465567), ("10", 97611390)]


# Test 1: vcf_get_sample_list
def test_vcf_get_sample_list_reads_header_samples_from_fixture(vcf_path):
    samples = vcf_get_sample_list(vcf_path)
    # Just sanity-check that we got some samples and the type is correct
    assert isinstance(samples, list)
    assert len(samples) > 0
    assert len(samples) == 3202

# Test 2: vcf_buffer_variants <- run_bcftools_view, parse_target_variant <- vcf_get_genotypes
# Test 2.1: normal case
def test_vcf_buffer_variants(test_positions, vcf_path, tmp_path_factory,sampleset_name):
    cache_dir = tmp_path_factory.mktemp("cache")
    variants = vcf_buffer_variants(
        position_batch=test_positions, 
        target_path=vcf_path, 
        cache_dir=cache_dir,
        sampleset=sampleset_name
        )
    
    # test output type
    assert isinstance(variants, TargetVariants)
    # test genotype properties
    assert variants.genotypes.shape == (3, 3202, 2)
    # test variant_metadata properties
    assert variants.variant_metadata.chr_name == ["1", "8", "10"]
    assert variants.variant_metadata.chr_pos == [44465567, 116938529, 97611390]
    assert variants.variant_metadata.ref == ['G', 'G', 'A']
    assert variants.variant_metadata.alts == [['A'], ['C'], ['G']]
    # test variant_ids
    assert variants.variant_ids == ['1:44465567:G:A', '8:116938529:G:C', '10:97611390:A:G']

# Test 2.2: missing positions
def test_vcf_buffer_variants_with_missing(test_positions, vcf_path, tmp_path_factory, sampleset_name, n_samples):
    missing_chrom, missing_pos = "1", 9067156  # definitely missing
    positions = test_positions + [(missing_chrom, missing_pos)]
    cache_dir = tmp_path_factory.mktemp("cache")
    variants = vcf_buffer_variants(
            position_batch=positions,
            target_path=vcf_path,
            cache_dir=cache_dir,
            sampleset=sampleset_name
        )
    
    assert len(variants) == len(positions)
    missing_index = positions.index((missing_chrom, missing_pos))
    assert variants.variant_metadata.chr_name[missing_index] == "1"
    assert variants.variant_metadata.ref[missing_index] is None
    assert variants.variant_metadata.alts[missing_index] is None

    # test that missing genotypes are filled with MISSING_GENOTYPE_TUPLE
    missing_gts = variants.genotypes[missing_index]
    print(missing_gts.shape)
    assert missing_gts.shape == (n_samples, 2)
    assert missing_gts.dtype == np.uint8
    assert np.all(missing_gts == MISSING_GENOTYPE_TUPLE)

# Test 2.3: vcf file without tabix index
def test_vcf_buffer_no_index(test_positions, vcf_path, tmp_path_factory, sampleset_name):
    """An unindexed VCF must throw a pysam.SamtoolsError"""
    vcf_dir = tmp_path_factory.mktemp("vcf")
    shutil.copy(vcf_path, vcf_dir)
    cache_dir = tmp_path_factory.mktemp("cache")

    with pytest.raises(pysam.SamtoolsError) as e:
        _ = vcf_buffer_variants(
                position_batch=test_positions,
                target_path=vcf_dir / vcf_path.name,
                cache_dir=cache_dir,
                sampleset=sampleset_name
            )
    assert "could not load index" in str(e.value)

# Test 2.4: with X/Y haploid genotypes
def test_vcf_buffer_sex_variants(x_y_variants, vcf_path, tmp_path_factory, sampleset_name):
    """Haploid genotypes will cause this function to throw NotImplementedErrors"""
    cache_dir = tmp_path_factory.mktemp("cache")
    with pytest.raises(pysam.SamtoolsError):
        _ = vcf_buffer_variants(
                position_batch=x_y_variants, 
                target_path=vcf_path, 
                cache_dir=cache_dir,
                sampleset=sampleset_name
            )