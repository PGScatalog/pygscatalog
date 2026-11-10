import dask.array as da
import numpy as np
import pytest
import zarr

from pgscatalog.calc.lib.constants import (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    NO_ALLELE_SENTINEL_VALUE,
)
from pgscatalog.calc.lib.score._dosage import (
    adjust_dosage_for_effect,
    calculate_effect_allele_dosage,
    fill_missing_dosages,
    no_called_alleles_reason,
    sentinel_to_nan,
)

N_SAMPLES = 5


@pytest.fixture
def dosage():
    # first variant: 0 copies of effect allele, second: 1 copy, third: 2 copies
    dosage = np.zeros(shape=(N_SAMPLES, 1000))
    dosage[1,] = 1
    dosage[2,] = 2
    dosage[3,] = 1
    dosage[4,] = 2
    return dosage


@pytest.fixture
def dosage_ones():
    return np.ones(shape=(N_SAMPLES, 1000))


def test_adjust_dosage_for_effect(dosage, tmp_path):
    is_recessive = np.array([False, True, True, False, False])
    is_dominant = np.array([False, False, False, True, True])
    store = zarr.storage.MemoryStore()
    dosage_zarr = zarr.create_array(store=store, shape=dosage.shape, dtype=dosage.dtype)
    dosage_zarr[:] = dosage
    adjusted = adjust_dosage_for_effect(
        dosage_array=dosage_zarr, recessive_mask=is_recessive, dominant_mask=is_dominant
    )

    # 0/0
    assert np.all(adjusted[0,] == 0)  # additive = no change

    # 0/1
    assert np.all(adjusted[1,] == 0)  # recessive = max(dosage - 1, 0)
    assert np.all(adjusted[3,] == 1)  # dominant = (greater than 1 is 1)

    # 1/1
    assert np.all(adjusted[2,] == 1)  # recessive = max(dosage - 1, 0)
    assert np.all(adjusted[4,] == 1)  # dominant = (greater than 1 is 1)


def test_bad_effect_types(dosage, tmp_path):
    # both recessive and dominant is impossible
    is_recessive = np.array([False, True, True, False, False])
    is_dominant = np.array([False, False, True, True, True])

    store = zarr.storage.MemoryStore()
    dosage_zarr = zarr.create_array(store=store, shape=dosage.shape, dtype=dosage.dtype)
    dosage_zarr[:] = dosage
    with pytest.raises(ValueError) as e:
        _ = adjust_dosage_for_effect(
            dosage_array=dosage_zarr,
            recessive_mask=is_recessive,
            dominant_mask=is_dominant,
        )
    assert "A variant cannot be both" in str(e.value)


@pytest.fixture
def db_metadata():
    return "PGS001229_hmPOS_GRCh38", "1000G"


def test_a_hemizygous_call_is_never_treated_as_missing():
    """The no-allele sentinel must never reach NaN (ADR 0007).

    `sentinel_to_nan` is the most dangerous line in haploid support: if 254 went
    down the missing path, `calculate_effect_allele_dosage` would NaN the whole
    call and every man's entire X chromosome would be imputed, silently, into a
    plausible number. Only 255 is a missing call.
    """
    hemi_alt, hemi_ref = 1, 0
    gts = np.array(
        [
            [
                [hemi_alt, NO_ALLELE_SENTINEL_VALUE],  # hemizygous ALT
                [hemi_ref, NO_ALLELE_SENTINEL_VALUE],  # hemizygous REF
                [0, 1],  # diploid het
                [1, 1],  # diploid hom ALT
                [MISSING_GENOTYPE_SENTINEL_VALUE] * 2,  # diploid missing
                [MISSING_GENOTYPE_SENTINEL_VALUE, NO_ALLELE_SENTINEL_VALUE],
            ]
        ],
        dtype=np.uint8,
    )
    masked = sentinel_to_nan(da.from_array(gts))

    # the sentinel survives the mask: nothing but 255 is a missing call
    assert not bool(da.isnan(masked[0, 0, 1]).compute())

    dosage = calculate_effect_allele_dosage(
        masked, effect_idx=np.array([1], dtype=np.int64)
    ).compute()

    # a hemizygous effect allele contributes 1, not 2, and is not missing
    np.testing.assert_array_equal(dosage[0], [1, 0, 1, 2, np.nan, np.nan])


def test_ploidy_is_derived_from_the_genotype_array():
    """Never stored: a mask would duplicate what the array already says (ADR 0007)"""
    gts = np.array(
        [
            [1, NO_ALLELE_SENTINEL_VALUE],
            [0, 1],
            [MISSING_GENOTYPE_SENTINEL_VALUE, NO_ALLELE_SENTINEL_VALUE],
            [MISSING_GENOTYPE_SENTINEL_VALUE] * 2,
        ],
        dtype=np.uint8,
    )

    ploidy = np.sum(gts != NO_ALLELE_SENTINEL_VALUE, axis=1)

    np.testing.assert_array_equal(ploidy, [1, 2, 1, 2])


def test_small_shard_without_missing_calls_skips_the_floor():
    """The sample floor only applies if there's something to impute, so a tiny
    shard with complete calls scores instead of crashing"""
    dosage = da.from_array(np.ones(shape=(4, 10)))
    ploidy = da.full(dosage.shape, 2, dtype=np.uint8)

    is_nan, filled = fill_missing_dosages(
        dosage_array=dosage,
        ploidy_array=ploidy,
        n_minimum_samples=50,
        effect_allele_frequencies=np.full(4, 0.5),
        n_samples=10,
    )

    assert not bool(is_nan.any().compute())
    np.testing.assert_array_equal(filled.compute(), dosage.compute())


def test_small_shard_with_missing_calls_still_raises():
    """...but imputing from 10 samples is still not allowed"""
    array = np.ones(shape=(4, 10))
    array[0, 0] = np.nan
    dosage = da.from_array(array)
    ploidy = da.full(dosage.shape, 2, dtype=np.uint8)

    with pytest.raises(ValueError) as e:
        fill_missing_dosages(
            dosage_array=dosage,
            ploidy_array=ploidy,
            n_minimum_samples=50,
            effect_allele_frequencies=np.full(4, 0.5),
            n_samples=10,
        )

    assert "n_minimum_samples=50" in str(e.value)


def test_a_missing_call_is_imputed_from_its_own_ploidy():
    """The imputed value is per (variant, sample), not a per-variant scalar: a
    man's missing X call is one expected copy, a woman's two (ADR 0009)"""
    #                  man  man  woman  woman  man(missing)  woman(missing)
    array = np.array([[1.0, 0.0, 2.0, 1.0, np.nan, np.nan]])
    ploidy = da.from_array(np.array([[1, 1, 2, 2, 1, 2]], dtype=np.uint8))

    _, filled = fill_missing_dosages(
        dosage_array=da.from_array(array),
        ploidy_array=ploidy,
        n_minimum_samples=1,
        effect_allele_frequencies=np.array([4 / 6]),
        n_samples=6,
    )

    # by hand: 4 effect alleles over 6 called ones
    freq = 4 / 6
    np.testing.assert_allclose(filled.compute()[0], [1, 0, 2, 1, freq, 2 * freq])


def _zarr(array):
    store = zarr.storage.MemoryStore()
    zarr_array = zarr.create_array(store=store, shape=array.shape, dtype=array.dtype)
    zarr_array[:] = array
    return zarr_array


def test_no_called_alleles_message_blames_missingness():
    """Everybody could carry it, but nobody was called"""
    dosage = _zarr(np.full((1, 4), np.nan))
    ploidy = _zarr(np.full((1, 4), 2, dtype=np.uint8))

    assert no_called_alleles_reason(dosage, ploidy) == "every call is missing"


def test_no_called_alleles_message_blames_nullploidy():
    """...and the other message, which nullploidy makes reachable (ADR 0007).

    Without it the first user to see this reads a biological fact as corruption
    """
    dosage = _zarr(np.full((1, 4), np.nan))
    ploidy = _zarr(np.zeros((1, 4), dtype=np.uint8))

    assert "nullploid" in no_called_alleles_reason(dosage, ploidy)
