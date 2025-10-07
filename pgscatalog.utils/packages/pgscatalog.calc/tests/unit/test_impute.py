import numpy as np
import pytest
import zarr

from pgscatalog.calc.lib._dosage import (
    fill_missing_dosages,
)
from pgscatalog.calc.lib._impute import calculate_mean_dosage
from pgscatalog.calc.lib.constants import MISSING_GENOTYPE_SENTINEL_VALUE


@pytest.fixture(scope="function")
def imputed_effect_dosage(missing_dosage_array):
    return calculate_mean_dosage(dosage_array=missing_dosage_array)


def test_good_impute(missing_dosage_array, imputed_effect_dosage):
    assert imputed_effect_dosage.shape[0] == missing_dosage_array.shape[0]
    assert np.all(0 <= imputed_effect_dosage) and np.all(2 >= imputed_effect_dosage), (
        "Dosage must be between 0 and 2"
    )


def test_fill_missing_with_imputed(missing_dosage_array, imputed_effect_dosage):
    nan_mask, filled = fill_missing_dosages(missing_dosage_array, n_minimum_samples=50)
    filled = filled.compute()

    # test that no NaNs remain
    assert not np.any(np.isnan(filled))

    # check that each NaN was replaced with the correct per-variant imputed value
    for variant_idx in range(missing_dosage_array.shape[0]):
        assert np.allclose(
            filled[variant_idx, nan_mask[variant_idx]],
            imputed_effect_dosage[variant_idx],
        )


def test_too_few_samples(missing_dosage_array):
    # imputation relies on estimating allelic frequency from non-missing samples
    # too few samples will return invalid dosages
    n_minimum_samples = missing_dosage_array.shape[1] + 1  # more than in dataset
    with pytest.raises(ValueError) as e:
        _ = fill_missing_dosages(
            dosage_array=missing_dosage_array, n_minimum_samples=n_minimum_samples
        )
    assert "is less than" in str(e.value)


def test_bad_type(missing_dosage_array):
    # missing values are represented using np.nan
    # integer types can't represent np.nan, so input must be floating
    store = zarr.storage.MemoryStore()
    bad_type = zarr.create_array(
        store=store,
        shape=missing_dosage_array.shape,
        dtype=np.uint8,
        fill_value=MISSING_GENOTYPE_SENTINEL_VALUE,
    )
    bad_type[:] = np.nan_to_num(
        missing_dosage_array, MISSING_GENOTYPE_SENTINEL_VALUE
    ).astype(np.uint8)
    with pytest.raises(TypeError) as e:
        _ = fill_missing_dosages(dosage_array=bad_type)
    assert "The genotype array must contain floats" in str(e.value)
