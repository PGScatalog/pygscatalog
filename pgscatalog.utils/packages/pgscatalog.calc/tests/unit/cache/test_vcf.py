import numpy as np
import pytest

from pgscatalog.calc.lib.cache._vcf import (
    region_sort_key,
    vcf_buffer_variants,
    vcf_get_sample_list,
)
from pgscatalog.calc.lib.cache.targetvariants import (
    TargetVariants,
)
from pgscatalog.calc.lib.constants import (
    MISSING_GENOTYPE_TUPLE,
    NO_ALLELE_SENTINEL_VALUE,
)


@pytest.fixture
def test_positions() -> list[tuple[str, int]]:
    """Some positions in tiny1000G"""
    return [("8", 116938529), ("1", 44465567), ("10", 97611390)]


def test_region_sort_key_accepts_non_numeric_chromosomes():
    """int("X") raised, so nothing downstream could run on a non-autosome"""
    regions = [("X", 2, 2), ("10", 1, 1), ("MT", 1, 1), ("2", 9, 9), ("2", 1, 1)]

    assert sorted(regions, key=region_sort_key) == [
        ("2", 1, 1),
        ("2", 9, 9),
        ("10", 1, 1),
        ("MT", 1, 1),
        ("X", 2, 2),
    ]


# Test 1: vcf_get_sample_list
def test_vcf_get_sample_list_reads_header_samples_from_fixture(vcf_path):
    samples = vcf_get_sample_list(vcf_path)
    # Just sanity-check that we got some samples and the type is correct
    assert isinstance(samples, list)
    assert len(samples) == 3202


# Test 2: vcf_buffer_variants <- run_bcftools_view,
#         parse_target_variant <- vcf_get_genotypes
# Test 2.1: normal case
def test_vcf_buffer_variants(
    test_positions, vcf_path, vcf_index_path, tmp_path_factory, sampleset_name
):
    cache_dir = tmp_path_factory.mktemp("cache")
    variants = vcf_buffer_variants(
        position_batch=test_positions,
        target_path=vcf_path,
        index_path=vcf_index_path,
        cache_dir=cache_dir,
        sampleset=sampleset_name,
    )

    # test output type
    assert isinstance(variants, TargetVariants)
    # test genotype properties
    assert variants.genotypes.shape == (3, 3202, 2)
    # test variant_metadata properties
    assert variants.variant_metadata.chr_name == ["1", "8", "10"]
    assert variants.variant_metadata.chr_pos == [44465567, 116938529, 97611390]
    assert variants.variant_metadata.ref == ["G", "G", "A"]
    assert variants.variant_metadata.alts == [["A"], ["C"], ["G"]]
    # test variant_ids
    assert variants.variant_ids == [
        "1:44465567:G:A",
        "8:116938529:G:C",
        "10:97611390:A:G",
    ]


# Test 2.2: missing positions
def test_vcf_buffer_variants_with_missing(
    test_positions,
    vcf_path,
    vcf_index_path,
    tmp_path_factory,
    sampleset_name,
    n_samples,
):
    missing_chrom, missing_pos = "1", 9067156  # definitely missing
    positions = test_positions + [(missing_chrom, missing_pos)]
    cache_dir = tmp_path_factory.mktemp("cache")
    variants = vcf_buffer_variants(
        position_batch=positions,
        target_path=vcf_path,
        index_path=vcf_index_path,
        cache_dir=cache_dir,
        sampleset=sampleset_name,
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


# Test 2.3: with haploid chrX genotypes
def test_vcf_buffer_x_variants(
    x_variants, vcf_path, vcf_index_path, tmp_path_factory, sampleset_name
):
    """Ploidy is read from the length of the GT field (ADR 0007).

    1000G writes a man's non-pseudoautosomal X call as `1`, not `1/1`, so the
    fixture holds 1599 hemizygous calls and 1603 diploid ones at every variant.
    Before this the array was built with `np.fromiter(..., dtype=(np.uint8, 2))`,
    which broadcasts `(1,)` to `[1, 1]` and quietly doubles every man's dosage.
    """
    cache_dir = tmp_path_factory.mktemp("cache")
    variants = vcf_buffer_variants(
        position_batch=x_variants,
        target_path=vcf_path,
        index_path=vcf_index_path,
        cache_dir=cache_dir,
        sampleset=sampleset_name,
    )

    assert set(variants.variant_metadata.chr_name) == {"X"}

    ploidy = np.sum(variants.genotypes != NO_ALLELE_SENTINEL_VALUE, axis=2)
    called = variants.genotypes[ploidy == 1]

    # a hemizygous call keeps its one allele: nothing was duplicated into the
    # slot the sample doesn't have
    assert np.all(np.isin(called[:, 0], [0, 1]))
    assert np.all(called[:, 1] == NO_ALLELE_SENTINEL_VALUE)

    # every variant the target actually holds is called in the same mixture of
    # men and women, which is why ploidy is a property of a call
    genotyped = [ref is not None for ref in variants.variant_metadata.ref]
    assert {tuple(np.bincount(row, minlength=3)[1:]) for row in ploidy[genotyped]} == {
        (1599, 1603)
    }
