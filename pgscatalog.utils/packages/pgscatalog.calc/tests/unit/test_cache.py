import json
import pathlib

import numpy as np
import pysam
import pytest
import zarr

from pgscatalog.calc import Scorefiles, TargetGenome
from pgscatalog.calc.lib.constants import (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    ZARR_MAX_N_VARIANTS,
    ZARR_PLOIDY,
)


@pytest.fixture(scope="session")
def vcf_cache(tmp_path_factory, pgs001229, vcf_path, sampleset_name) -> pathlib.Path:
    tmpdir = tmp_path_factory.mktemp("vcf_cache")
    genome = TargetGenome(vcf_path, cache_dir=tmpdir, sampleset=sampleset_name)
    scorefile = Scorefiles(pgs001229)
    positions = scorefile.get_unique_positions()
    genome.cache_variants(positions=positions)
    return pathlib.Path(tmpdir)


def test_cache_structure(vcf_cache, sampleset_name, vcf_path):
    """Test the process of turning a VCF file into a cache"""
    # test hive partitioning structure of zarr files
    assert {x.name for x in (vcf_cache / "gts").iterdir()} == {"zarr.json", "test"}
    assert {x.name for x in (vcf_cache / "gts" / sampleset_name).iterdir()} == {
        "zarr.json",
        vcf_path.name,
    }
    assert (vcf_cache / "variants.db").exists()


def test_cache_genotypes(vcf_cache, n_samples, sampleset_name, vcf_path):
    """Test the genotype array in a zarr local store"""
    store = zarr.storage.LocalStore(vcf_cache / "gts")
    root = zarr.open_group(store=store, path=sampleset_name)
    array = root[vcf_path.name]["genotypes"]

    assert array.shape == (ZARR_MAX_N_VARIANTS, n_samples, ZARR_PLOIDY)
    assert array.dtype == np.uint8

    # simple checks for the genotype array
    # the placeholder values have been replaced (zarr array is initialised with np.nan)
    assert not np.all(np.isnan(array[:1000, :, :]))
    # check that not every genotype is missing
    assert not np.all(array[:1000, :, :] == MISSING_GENOTYPE_SENTINEL_VALUE)


def test_cache_metadata(vcf_cache, vcf_path, sampleset_name):
    sampleset_metadata = vcf_cache / "gts" / sampleset_name / "zarr.json"

    with (
        sampleset_metadata.open() as f,
        pysam.VariantFile(str(vcf_path), mode="r") as v,
    ):
        metadata = json.load(f)
        assert metadata["attributes"]["samples"] == list(v.header.samples)
