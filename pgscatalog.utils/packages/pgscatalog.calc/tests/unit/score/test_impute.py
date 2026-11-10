import numpy as np
import pytest
from dask import array as da

from pgscatalog.calc.lib.score._dosage import fill_missing_dosages


def test_small_shard_uses_sampleset_frequency_and_sample_floor():
    # No called alleles here: all evidence came from the other shards.
    dosage = da.from_array(np.array([[np.nan, np.nan, 0.0]]))
    ploidy = da.from_array(np.array([[1, 2, 0]], dtype=np.uint8))
    missing, filled = fill_missing_dosages(
        dosage,
        ploidy,
        np.array([0.25]),
        n_samples=50,
    )
    np.testing.assert_array_equal(missing.compute(), [[True, True, False]])
    np.testing.assert_array_equal(filled.compute(), [[0.25, 0.5, 0]])


def test_all_nullploid_calls_need_no_frequency_or_sample_floor():
    dosage = da.zeros((1, 3))
    _, filled = fill_missing_dosages(
        dosage,
        da.zeros((1, 3), dtype=np.uint8),
        np.array([np.nan]),
        n_samples=3,
    )
    np.testing.assert_array_equal(filled.compute(), [[0, 0, 0]])


def test_sampleset_below_floor_cannot_impute():
    with pytest.raises(ValueError, match="n_samples=49"):
        fill_missing_dosages(
            da.from_array(np.array([[np.nan]])),
            da.ones((1, 1), dtype=np.uint8),
            np.array([0.5]),
            n_samples=49,
        )
