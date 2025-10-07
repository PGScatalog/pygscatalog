import numpy as np
import numpy.typing as npt
import pytest

from pgscatalog.calc import TargetVariant, TargetVariants


@pytest.fixture
def variant() -> TargetVariant:
    gts: npt.NDArray[np.uint8] = np.array([[0, 1], [1, 1], [0, 0]], dtype=np.uint8)
    return TargetVariant(chr_name="21", chr_pos=46366609, ref="A", alts=("G",), gts=gts)


@pytest.fixture
def padded_variant() -> TargetVariant:
    return TargetVariant(chr_name="01", chr_pos=46366609, ref=None, alts=None, gts=None)


def test_variant():
    variant = TargetVariant(chr_name="1", chr_pos=1, ref=None, alts=None, gts=None)
    assert variant.chr_name == "1"


def test_padded_variant(padded_variant):
    assert padded_variant.chr_name == "1"


def test_bad_variant(variant) -> None:
    # test wrong number of sample IDs fails properly
    with pytest.raises(ValueError) as excinfo:
        TargetVariants(
            variants=[variant],
            samples=["test"],
            filename="test.vcf.gz",
            sampleset="test",
        )
    assert "Number of samples in genotypes" in str(excinfo.value)

    # test duplicated sample IDs fail properly
    with pytest.raises(ValueError) as excinfo:
        TargetVariants(
            variants=[variant],
            samples=["test", "test", "test"],
            filename="test.vcf.gz",
            sampleset="test",
        )
    assert "Number of samples in genotypes" in str(excinfo.value)

    # badly shaped genotype array
    gts: npt.NDArray[np.uint8] = np.array([[0, 1, 1], [0, 0, 0]], dtype=np.uint8)
    with pytest.raises(ValueError) as excinfo:
        TargetVariant(chr_name="21", chr_pos=46366609, ref="A", alts=("G",), gts=gts)
    assert "Array shape" in str(excinfo.value)
