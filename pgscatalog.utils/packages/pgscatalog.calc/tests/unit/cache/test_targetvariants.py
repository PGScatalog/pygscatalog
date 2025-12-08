import numpy as np
import pytest
import zarr

from pgscatalog.calc.lib.cache.targetvariants import (
    TargetVariants,
    add_missing_positions_to_lists,
)
from pgscatalog.calc.lib.cache.zarrmodels import ZarrSampleMetadata, ZarrVariantMetadata
from pgscatalog.calc.lib.constants import MISSING_GENOTYPE_SENTINEL_VALUE


@pytest.fixture
def valid_args(tmp_path) -> dict:
    """Common valid constructor args for TargetVariants."""
    chr_name = ["1", "21"]
    chr_pos = [46366609, 46366609]
    refs = ["A", "C"]
    alts = [["G"], ["A"]]

    gts = [
        np.array([[0, 1], [1, 1], [0, 0]], dtype=np.uint8),
        np.array(
            [
                [0, 0],
                [0, 1],
                [MISSING_GENOTYPE_SENTINEL_VALUE, MISSING_GENOTYPE_SENTINEL_VALUE],
            ],
            dtype=np.uint8,
        ),
    ]
    samples = ["sample1", "sample2", "sample3"]
    target_path = tmp_path / "dummy.bgen"

    return {
        "chr_name": chr_name,
        "pos": chr_pos,
        "refs": refs,
        "alts": alts,
        "gts": gts,
        "samples": samples,
        "sampleset": "valid",
        "target_path": target_path,
    }


@pytest.fixture
def variant(valid_args) -> TargetVariants:
    """A valid TargetVariants instance for most tests."""
    return TargetVariants(**valid_args)


# Test 1: Test properties of the minimal construction
def test_variant_basic_construction(variant: TargetVariants):
    test_variant = variant
    # test genotype properties
    assert test_variant.genotypes.shape == (2, 3, 2)
    assert test_variant.genotypes.dtype == np.uint8

    # test samples property
    assert test_variant.samples == ["sample1", "sample2", "sample3"]

    # test variant IDs property
    ids = test_variant.variant_ids
    assert len(ids) == len(test_variant)
    assert ids[0] == "1:46366609:A:G"

    # test variant_metadata property
    sample_metadata = test_variant.variant_metadata
    assert isinstance(sample_metadata, ZarrVariantMetadata)
    assert sample_metadata.chr_name[0] == "1"
    assert sample_metadata.chr_pos[0] == 46366609
    assert sample_metadata.ref[0] == "A"
    assert sample_metadata.alts[0] == ["G"]


# Test 2: Test that bad constructions raise errors
# error 1: mismatched_lengths among inputs
def test_variant_bad_construction(valid_args: dict):
    mismatched_args = valid_args.copy()
    # drop one alt so lengths differ
    mismatched_args["alts"] = [["G"]]

    with pytest.raises(ValueError, match="different lengths"):
        TargetVariants(**mismatched_args)


# error 2: invalid values in genotypes
def test_variant_invalid_genotypes_construction(valid_args: dict):
    invalid_gts_args = valid_args.copy()
    # introduce an invalid genotype value 2
    invalid_gts = [
        valid_args["gts"][0],
        np.array([[0, 0], [0, 1], [2, 1]], dtype=np.uint8),  # invalid '2'
    ]
    invalid_gts_args["gts"] = invalid_gts

    with pytest.raises(ValueError, match="Genotypes contain invalid values"):
        TargetVariants(**invalid_gts_args)


# Test 3: test_write_zarr_writes_expected_structure
def _make_memory_group() -> zarr.Group:
    """Create an in-memory zarr root group for testing."""
    root = zarr.group()
    return root


# Zarr 1. build the zarr
def test_write_zarr(variant: TargetVariants):
    root = _make_memory_group()
    variant.write_zarr(root)

    # 1) genotypes array
    assert "genotypes" in root
    gt_arr = root["genotypes"]
    assert gt_arr.shape == variant.genotypes.shape
    assert gt_arr.dtype == np.uint8
    np.testing.assert_array_equal(gt_arr[...], variant.genotypes)

    # 2) samples attribute
    assert "samples" in root.attrs
    assert (
        root.attrs["samples"] == ZarrSampleMetadata.model_validate(variant.samples).root
    )

    # 3) meta subgroup with 1D arrays
    meta_root = zarr.open_group(store=root.store, path="meta")
    for key in ("chr_name", "chr_pos", "ref", "alts", "variant_id"):
        assert key in meta_root


# Zarr 2. append when arrays exist
def test_write_zarr_appends_when_arrays_exist(variant: TargetVariants):
    root = _make_memory_group()

    # first write
    variant.write_zarr(root)
    n_variants = len(variant)

    # second write should append rows
    variant.write_zarr(root)

    # check genotypes shape doubled
    gt_arr = root["genotypes"]
    assert gt_arr.shape[0] == 2 * n_variants

    # check meta arrays variant_id doubled
    meta = zarr.open_group(store=root.store, path="meta")
    variant_id_arr = meta["variant_id"][...]
    variant_id_values = [
        x.decode() if isinstance(x, bytes) else x for x in variant_id_arr
    ]
    assert variant_id_values == variant.variant_ids * 2


# Test 4: Test add_missing_positions_to_lists function
def test_add_missing_positions_to_lists(variant: TargetVariants):
    chroms = ["1"]
    positions = [46366609]
    ref_alleles = ["A"]
    alt_alleles = [["G"]]

    existing_gts = np.array([[0, 1], [1, 0]], dtype=np.uint8)
    hard_calls = [existing_gts]

    # ("2", 46366610) is missing
    existing_positions = {("1", 46366609)}
    scoring_file_regions = [("1", 46366609), ("2", 46366610)]
    n_samples = 2

    add_missing_positions_to_lists(
        chroms=chroms,
        positions=positions,
        ref_alleles=ref_alleles,
        alt_alleles=alt_alleles,
        hard_calls=hard_calls,
        scoring_file_regions=scoring_file_regions,
        seen_positions=existing_positions,
        n_samples=n_samples,
    )

    # one extra position added
    assert len(chroms) == 2
    assert len(positions) == 2
    assert len(ref_alleles) == 2
    assert len(alt_alleles) == 2
    assert len(hard_calls) == 2

    # find missing position
    missing_index = positions.index(46366610)
    assert chroms[missing_index] == "2"
    assert ref_alleles[missing_index] is None
    assert alt_alleles[missing_index] is None

    missing_gts = hard_calls[missing_index]
    assert missing_gts.shape == (n_samples, 2)
    assert missing_gts.dtype == np.uint8
    assert np.all(missing_gts == MISSING_GENOTYPE_SENTINEL_VALUE)
