# import numpy as np
# import pytest
# import zarr
# from pgscatalog.calc.lib._bgen_probabilities import (
#     phased_probabilities_to_hard_calls,
#     unphased_probabilities_to_hard_calls,
# )
#
# from pgscatalog.calc import TargetGenome
# from pgscatalog.calc.lib.constants import MISSING_GENOTYPE_SENTINEL_VALUE
#
#
# @pytest.mark.parametrize(
#     ("target_with_missingness", "missing_bgen_sample"),
#     [
#         ("target_with_missingness", None),
#         ("missing_bgen", "missing_bgen_sample"),
#     ],
#     indirect=True,
#     ids=["VCF missingness", "bgen unphased missingness"],
# )
# def test_missing_encoded_correctly(
#     target_with_missingness,
#     test_positions,
#     tmp_path_factory,
#     missing_prob,
#     missing_bgen_sample,
#     sampleset_name,
# ) -> None:
#     cache_dir = tmp_path_factory.mktemp("cache")
#     TargetGenome(
#         target_path=target_with_missingness,
#         cache_dir=cache_dir,
#         sampleset=sampleset_name,
#     ).cache_variants(positions=test_positions)
#
#     store = zarr.storage.LocalStore(cache_dir / "gts")
#     group = zarr.open_group(
#         store=store, path=f"{sampleset_name}/{target_with_missingness.name}"
#     )
#     array = group["genotypes"]
#
#     gt_array = array[: len(test_positions), :, :]
#     n_missing = np.sum(gt_array == MISSING_GENOTYPE_SENTINEL_VALUE)
#
#     # a bit of wiggle room
#     min_missing = missing_prob - 0.02
#     assert float(n_missing / gt_array.size) > min_missing
#
#
# @pytest.fixture
# def phased_missing(missing_prob):
#     np.random.seed(42)
#     n_rows = 5000
#     n_cols = 2
#
#     # apparently a dirichlet distribution sums to 1
#     hap1 = np.random.dirichlet(np.ones(n_cols), size=n_rows)
#     hap2 = np.random.dirichlet(np.ones(n_cols), size=n_rows)
#
#     n_nan_rows = int(n_rows * missing_prob)
#     nan_row_indices = np.random.choice(n_rows, size=n_nan_rows, replace=False)
#
#     hap1[nan_row_indices] = np.nan
#     hap2[nan_row_indices] = np.nan
#
#     phased = []
#     for x in np.hstack([hap1, hap2]):
#         phased.append(np.atleast_2d(x))
#     return phased
#
#
# @pytest.fixture
# def unphased_missing(missing_prob):
#     np.random.seed(42)
#     n_rows = 5000
#     n_cols = 3
#
#     # apparently a dirichlet distribution sums to 1
#     probs = np.random.dirichlet(np.ones(n_cols), size=n_rows)
#
#     n_nan_rows = int(n_rows * missing_prob)
#     nan_row_indices = np.random.choice(n_rows, size=n_nan_rows, replace=False)
#
#     probs[nan_row_indices] = np.nan
#
#     unphased = []
#     for x in probs:
#         unphased.append(np.atleast_2d(x))
#
#     return unphased
#
#
# def test_phased_probabilities_to_hard_calls(phased_missing):
#     """Test that np.nan probabilities are correctly handled"""
#     arr = phased_probabilities_to_hard_calls(phased_missing).compute()
#     values, counts = np.unique(arr, return_counts=True)
#
#     assert arr.max() == 255
#     assert np.array_equal(values, np.array([0, 1, 255], dtype=np.uint8))
#     assert np.array_equal(counts, np.array([3483, 3517, 3000]))
#
#
# def test_unphased_probabilities_to_hard_calls(unphased_missing):
#     arr = unphased_probabilities_to_hard_calls(unphased_missing).compute()
#     values, _ = np.unique(arr, return_counts=True)
#
#     assert arr.max() == 255
#     assert np.array_equal(values, np.array([0, 1, 255], dtype=np.uint8))
#
#
# def test_failing_impute_missing(gt_array_with_missingness):
#     # TODO: imputation should fail with too few samples
#     pass
#
#
# def test_impute_missing(
#     gt_array_with_missingness,
# ):
#     # TODO: imputation should work
#     pass
