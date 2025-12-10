import numpy as np
import pytest
import zarr

from pgscatalog.calc.lib.score._dosage import (
    adjust_dosage_for_effect,
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
