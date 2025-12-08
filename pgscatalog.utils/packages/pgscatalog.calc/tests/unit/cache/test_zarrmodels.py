import pytest
import polars as pl
import numpy as np

from pgscatalog.calc.lib.cache.zarrmodels import ZarrVariantMetadata, ZarrSampleMetadata
from pgscatalog.calc.lib.constants import NUMPY_STRING_DTYPE


# Test 1: invalid allele, validator should raise ValueError
def test_variant_metadata_invalid_ref_allele_raises():
    with pytest.raises(ValueError, match="Invalid allele"):
        ZarrVariantMetadata(
            chr_name=["1"],
            chr_pos=[1],
            ref=["N"],  # invalid allele
            alts=["A"],
            variant_id=["v1"],
        )


# Test 2: different lengths of input lists
def test_variant_metadata_even_length_error():
    with pytest.raises(ValueError, match="All list fields must have the same length"):
        ZarrVariantMetadata(
            chr_name=["1", "1", "X"],
            chr_pos=[1, 2],  # shorter
            ref=["A", "T", "G"],
            alts=[["C"], ["TG"], ["T"]],
            variant_id=["var1", "var2", "var3"],
        )


# Test 3: returned type check for ZarrVariantMetadata
def simple_test_data():
    return ZarrVariantMetadata(
        chr_name=["1", "1", "22"],
        chr_pos=[1, 2, 3],
        ref=["A", "T", "G"],
        alts=[["C"], ["TG"], ["T"]],
        variant_id=["var1", "var2", "var3"],
    )


def test_zarr_variant_metadata_type():
    zarr_variant_metadata = simple_test_data()
    assert isinstance(zarr_variant_metadata, ZarrVariantMetadata)
    assert zarr_variant_metadata.__len__() == 3


# Test 4: merge two ZarrVariantMetadata instances
def test_zarr_variant_metadata_merge():
    zarr_variant_metadata_1 = simple_test_data()
    zarr_variant_metadata_2 = ZarrVariantMetadata(
        chr_name=["X"],
        chr_pos=[4],
        ref=["G"],
        alts=[["A"]],
        variant_id=["var4"],
    )
    merged = zarr_variant_metadata_1.merge(zarr_variant_metadata_2)
    assert merged.__len__() == 4
    assert merged.chr_name == ["1", "1", "22", "X"]
    assert merged.chr_pos == [1, 2, 3, 4]
    assert merged.ref == ["A", "T", "G", "G"]
    assert merged.alts == [["C"], ["TG"], ["T"], ["A"]]
    assert merged.variant_id == ["var1", "var2", "var3", "var4"]

    # originals should not be mutated
    assert zarr_variant_metadata_2.chr_name == ["X"]


# Test 5: ZarrVariantMetadata into df
def test_variant_metadata_to_df():
    zarr_variant_metadata = simple_test_data()
    df = zarr_variant_metadata.to_df()

    assert isinstance(df, pl.DataFrame)
    assert df.shape == (3, 5)
    assert df.columns == ["chr_name", "chr_pos", "ref", "alts", "variant_id"]


# Test 6: ZarrVariantMetadata into numpy arrays
def test_variant_metadata_to_numpy():
    zarr_variant_metadata = simple_test_data()
    numpy_data = zarr_variant_metadata.to_numpy()

    assert isinstance(numpy_data, dict)
    assert set(numpy_data.keys()) == {
        "chr_name",
        "chr_pos",
        "ref",
        "alts",
        "variant_id",
    }
    for k, v in numpy_data.items():
        assert isinstance(v, np.ndarray)
        assert len(v) == len(zarr_variant_metadata)
        if k == "chr_pos":
            assert numpy_data[k].dtype == np.int64
        else:
            assert numpy_data[k].dtype == NUMPY_STRING_DTYPE


# Test 7: Test ZarrSampleMetadata length
def test_sample_metadata_len():
    sample_meta = ZarrSampleMetadata(root=["sample1", "sample2"])
    assert len(sample_meta) == 2
    assert sample_meta.root == ["sample1", "sample2"]
