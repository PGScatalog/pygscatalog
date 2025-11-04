import shutil

import pysam
import pytest

from pgscatalog.calc.lib._vcf import vcf_buffer_variants, vcf_get_sample_list


def test_vcf_buffer_variants(pgs001229_autosome, vcf_path, tmp_path_factory):
    """Test that variants can be extracted from VCF files"""
    cache_dir = tmp_path_factory.mktemp("cache")
    variants = list(
        vcf_buffer_variants(
            position_batch=pgs001229_autosome, target_path=vcf_path, cache_dir=cache_dir
        )
    )

    n_samples: set[int] = {len(vcf_get_sample_list(vcf_path))}
    buffer_n_samples: set[int] = {x.gts.shape[0] for x in variants if x.gts is not None}

    assert n_samples == buffer_n_samples
    assert {x.chr_name for x in variants} == {str(x) for x in list(range(1, 23))}


def test_vcf_with_missing(test_positions, vcf_path, tmp_path_factory):
    missing_chrom, missing_pos = "1", 9067156  # definitely missing
    positions = test_positions + [(missing_chrom, missing_pos)]
    cache_dir = tmp_path_factory.mktemp("cache")

    variants = list(
        vcf_buffer_variants(
            position_batch=positions,
            target_path=vcf_path,
            cache_dir=cache_dir,
        )
    )

    # test that missing variants are yielded properly
    assert len(variants) == len(positions)
    missing_variant = variants[-1]
    assert missing_variant.gts is None
    assert missing_variant.ref is None
    assert missing_variant.alts is None
    assert missing_variant.chr_name == missing_chrom
    assert missing_variant.chr_pos == missing_pos


def test_vcf_buffer_sex_variants(x_y_variants, vcf_path, tmp_path_factory):
    """Haploid genotypes will cause this function to throw NotImplementedErrors"""
    cache_dir = tmp_path_factory.mktemp("cache")
    with pytest.raises(NotImplementedError):
        _ = list(
            vcf_buffer_variants(
                position_batch=x_y_variants, target_path=vcf_path, cache_dir=cache_dir
            )
        )


def test_vcf_buffer_no_index(pgs001229_autosome, vcf_path, tmp_path_factory):
    """An unindexed VCF must throw a pysam.SamtoolsError"""
    vcf_dir = tmp_path_factory.mktemp("vcf")
    shutil.copy(vcf_path, vcf_dir)
    cache_dir = tmp_path_factory.mktemp("cache")

    with pytest.raises(pysam.SamtoolsError) as e:
        _ = list(
            vcf_buffer_variants(
                position_batch=pgs001229_autosome,
                target_path=vcf_dir / vcf_path.name,
                cache_dir=cache_dir,
            )
        )
    assert "could not load index" in str(e.value)
