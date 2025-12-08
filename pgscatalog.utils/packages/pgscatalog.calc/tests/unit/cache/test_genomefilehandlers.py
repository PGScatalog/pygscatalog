import pathlib
import pytest

from pgscatalog.calc.lib.cache._genomefilehandlers import (
    get_file_handler, 
    GenomeFileHandler, 
    VCFHandler, 
    BgenFileHandler
    )

from pgscatalog.calc.lib.cache.genomefiletypes import GenomeFileType
from pgscatalog.calc.lib.cache.targetvariants import TargetVariants
# Test 1 – VCF path ends with .vcf.gz → returns VCFHandler
def test_get_file_handler_vcf(vcf_path, tmp_path):

    handler = get_file_handler(
        path=vcf_path,
        cache_dir=tmp_path,
        sampleset="TEST",
    )

    assert isinstance(handler, VCFHandler)
    assert handler.file_type == GenomeFileType.VCF
    assert handler.target_path.name == "tiny1000G.vcf.gz"
    assert handler.index_path.name == "tiny1000G.vcf.gz.tbi"
    assert len(handler.chroms) > 0
    assert len(handler.samples) > 0

def test_vcf_handler_query_variants(
    vcf_path, 
    tmp_path, 
    test_positions):

    handler = get_file_handler(
        path=vcf_path,
        cache_dir=tmp_path,
        sampleset="TEST",
    )

    variants = handler.query_variants(test_positions)

    # basic sanity checks on the TargetVariants result
    assert isinstance(variants, TargetVariants)
    assert len(variants.variant_metadata.chr_name) == len(test_positions)
    assert len(variants.variant_metadata.chr_pos) == len(test_positions)
    # at least one genotype array present
    assert variants.genotypes is not None
    assert len(variants.genotypes) > 0

# Test 2.1: bgen + sample_file + index → BgenFileHandler
def test_get_file_handler_bgen_with_sample_file(
    phased_bgen_path, bgen_sample, tmp_path):

    handler = get_file_handler(
        path=phased_bgen_path,
        cache_dir=tmp_path,
        sampleset="TEST",
        sample_file=bgen_sample,
    )

    assert isinstance(handler, BgenFileHandler)
    assert handler.file_type == GenomeFileType.BGEN
    assert handler.target_path.name == "tiny1000G.bgen"
    assert handler.index_path.name.endswith(".bgen.bgi")
    assert len(handler.chroms) > 0
    assert len(handler.samples) > 0

def test_bgen_handler_query_variants(
    phased_bgen_path, 
    bgen_sample, 
    tmp_path, 
    test_positions):

    handler = get_file_handler(
        path=phased_bgen_path,
        cache_dir=tmp_path,
        sampleset="TEST",
        sample_file=bgen_sample,
    )

    # Before query, index should be the original .bgi
    original_index = handler._index_path

    variants = handler.query_variants(test_positions)

    # After query, index should have been copied into cache/tmp
    cache_index = tmp_path / "tmp" / original_index.name
    assert cache_index.exists()
    assert handler._index_path == cache_index

    # Basic sanity check on returned variants
    assert isinstance(variants, TargetVariants)
    assert len(variants.variant_metadata.chr_name) == len(test_positions)
    assert len(variants.variant_metadata.chr_pos) == len(test_positions)

# Test 2.2: bgen is provided. sample_file does not exist, it should raise a FileNot
def test_get_file_handler_bgen_without_sample_file(unphased_bgen_path, tmp_path):

    with pytest.raises(ValueError, match="BGEN files require a sample file"):
        get_file_handler(
            path=unphased_bgen_path,
            cache_dir=tmp_path,
            sampleset="TEST",
            sample_file=None, 
        )

# Test 3 – path does not exist → FileNotFoundError
def test_get_file_handler_file_not_found(tmp_path):
    # the input file does not exist
    missing_vcf = tmp_path / "i_do_not_exist.vcf.gz"

    with pytest.raises(FileNotFoundError):
        get_file_handler(
            path=missing_vcf,
            cache_dir=tmp_path,
            sampleset="TEST",
        )
    
    missing_bgen = tmp_path / "i_do_not_exist.bgen"
    missing_bgen_sample_file = tmp_path / "i_do_not_exist.sample"

    with pytest.raises(FileNotFoundError):
        get_file_handler(
            path=missing_bgen,
            cache_dir=tmp_path,
            sample_file=missing_bgen_sample_file,
            sampleset="TEST",
        )
    
# Test 4 – unsupported extension
def test_get_file_handler_unsupported_extension(tmp_path):
    fake_input = tmp_path / "weird.txt"
    fake_input.write_text("fake")

    with pytest.raises(ValueError, match="Unsupported genome file type"):
        get_file_handler(
            path=fake_input,
            cache_dir=tmp_path,
            sampleset="TEST",
        )

def test_bgen_copy_index_to_cache_creates_and_reuses_cached_index(
    phased_bgen_path,
    bgen_sample,
    tmp_path,
    caplog,
):
    handler = get_file_handler(
        path=phased_bgen_path,
        cache_dir=tmp_path,
        sampleset="TEST",
        sample_file=bgen_sample,
    )

    original_index = handler._index_path
    cache_index_path = tmp_path / "tmp" / original_index.name

    # First call: file does not exist yet → copy
    assert not cache_index_path.exists()

    with caplog.at_level("INFO"):
        handler.copy_index_to_cache()

    assert cache_index_path.exists()
    assert handler._index_path == cache_index_path

    # Second call: file exists → reuse branch
    caplog.clear()
    with caplog.at_level("INFO"):
        handler.copy_index_to_cache()

    # Ensure we hit the "already in cache" log message
    assert any("bgi index already in cache" in msg for msg in caplog.messages)
    assert handler._index_path == cache_index_path