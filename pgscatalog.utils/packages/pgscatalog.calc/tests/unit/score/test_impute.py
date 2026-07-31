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
    ploidy = da.full_like(dask_array, 2, dtype=np.uint8)

    result = calculate_mean_dosage(
        dask_array, ploidy_array=ploidy, n_minimum_samples=1
    ).compute()  # [0.375,0.5]
    # what we mathematically expect if missing samples are ignored
    expected = np.array(
        [
            ((0.0 + 1.0 + 2.0 + 0.0) / (3 * 2)) * 2,
            ((2.0 + 2.0 + 0.0 + 0.0) / (2 * 2)) * 2,
        ],
        dtype=np.float64,
    )

    assert result.shape == base.shape
    # every call is diploid here, so every sample gets the same imputed dosage
    assert np.allclose(result, expected[:, np.newaxis])


def test_calculate_mean_dosage_accepts_zarr_array(missing_dosage_array):
    """
    Check that if the calculate_mean_dosage works when passed the zarr-backed
    dosage array produced in the pipeline (missing_dosage_array fixture).
    ~30% of samples have randomly missing genotypes in each variant
    """
    ploidy = da.full(missing_dosage_array.shape, 2, dtype=np.uint8)
    result = calculate_mean_dosage(
        missing_dosage_array, ploidy_array=ploidy, n_minimum_samples=1
    )
    result_np = result.compute()

    assert result_np.shape == missing_dosage_array.shape
    # just a soft sanity check: values should be finite and within range
    assert np.all(np.isfinite(result_np))
    assert np.all(result_np >= 0)
    assert np.all(result_np <= 2)


def test_autosomal_imputation_is_bit_identical(missing_dosage_array):
    """Nothing moves on an all-diploid shard, to the last bit.

    `old` is the pre-0021 arithmetic written out: a per-variant scalar from a
    ploidy-blind denominator, broadcast across samples.
    """
    dosage = da.from_zarr(missing_dosage_array)
    ploidy = da.full(dosage.shape, 2, dtype=np.uint8)

    n_effect_alleles = da.nansum(dosage, axis=1)
    n_called_alleles = da.sum(~da.isnan(dosage), axis=1) * 2
    old = ((n_effect_alleles / n_called_alleles) * 2).compute()

    new = calculate_mean_dosage(
        dosage, ploidy_array=ploidy, n_minimum_samples=1
    ).compute()

    assert np.all(new == old[:, np.newaxis])


def test_the_frequency_is_pooled_across_ploidy():
    """A chrX variant's allele frequency is not halved by hemizygous samples.

    The expected frequency is computed by hand over the same calls, not from the
    implementation: 4 effect alleles over 7 called ones (ADR 0009).
    """
    #        man  man  woman  woman  woman(missing)
    dosage = np.array([[1.0, 0.0, 2.0, 1.0, np.nan]])
    ploidy = np.array([[1, 1, 2, 2, 2]], dtype=np.uint8)

    # by hand: effect alleles 1 + 0 + 2 + 1 = 4, called alleles 1 + 1 + 2 + 2 = 6
    freq = 4 / 6

    result = calculate_mean_dosage(
        da.from_array(dosage),
        ploidy_array=da.from_array(ploidy),
        n_minimum_samples=1,
    ).compute()

    # a ploidy-blind denominator would have used 4 called genotypes * 2 = 8 alleles,
    # so the frequency would have come out at 0.5 instead
    assert freq != 0.5
    np.testing.assert_allclose(result[0], [freq, freq, 2 * freq, 2 * freq, 2 * freq])


def test_a_nullploid_sample_imputes_to_zero_and_contributes_no_alleles():
    """A woman is never assigned a chrY genotype (ADR 0009)"""
    dosage = np.array([[1.0, 1.0, np.nan]])
    ploidy = np.array([[1, 1, 0]], dtype=np.uint8)

    result = calculate_mean_dosage(
        da.from_array(dosage),
        ploidy_array=da.from_array(ploidy),
        n_minimum_samples=1,
    ).compute()

    # both called alleles carry the effect allele, so the frequency is 1: the
    # nullploid sample is out of the denominator as well as out of the numerator
    np.testing.assert_array_equal(result[0], [1.0, 1.0, 0.0])
