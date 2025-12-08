import dask.array as da
import numpy as np

from pgscatalog.calc.lib.score._impute import calculate_mean_dosage


def test_calculate_mean_dosage_expected_values():
    base = np.array(
        [
            [0.0, 1.0, 2.0, np.nan],  # mean after drop na: (0 + 1 + 2) / 3*2 = 0.5
            [2.0, 2.0, np.nan, np.nan],  # mean after drop na: (2 + 2) / 2*2 = 1.0
        ],
        dtype=np.float64,
    )
    dask_array = da.from_array(base, chunks=(2, 2))

    result = calculate_mean_dosage(
        dask_array, n_minimum_samples=1
    ).compute()  # [0.375,0.5]
    # what we mathematically expect if missing samples are ignored
    expected = np.array(
        [
            ((0.0 + 1.0 + 2.0 + 0.0) / (3 * 2)) * 2,
            ((2.0 + 2.0 + 0.0 + 0.0) / (2 * 2)) * 2,
        ],
        dtype=np.float64,
    )

    assert result.shape == (base.shape[0],)
    assert np.allclose(result, expected)


def test_calculate_mean_dosage_accepts_zarr_array(missing_dosage_array):
    """
    Check that if the calculate_mean_dosage works when passed the zarr-backed
    dosage array produced in the pipeline (missing_dosage_array fixture).
    ~30% of samples have randomly missing genotypes in each variant
    """
    result = calculate_mean_dosage(missing_dosage_array, n_minimum_samples=1)
    result_np = result.compute()

    assert result_np.ndim == 1
    assert result_np.shape[0] == missing_dosage_array.shape[0]
    # just a soft sanity check: values should be finite and within range
    assert np.all(np.isfinite(result_np))
    assert np.all(result_np >= 0)
    assert np.all(result_np <= 2)
