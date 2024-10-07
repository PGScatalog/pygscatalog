import pytest
import csv
from pgscatalog.core.lib.models import CatalogScoreVariant


@pytest.fixture(scope="package")
def complex_variants(request):
    pgs000021 = request.path.parent / "data" / "complex_variants_PGS000021.txt"
    pgs004239 = request.path.parent / "data" / "complex_variants_PGS004239.txt"
    with open(pgs000021, "rt") as f1, open(pgs004239, "rt") as f2:
        return list(csv.DictReader(f1, delimiter="\t")), list(
            csv.DictReader(f2, delimiter="\t")
        )


def test_complex_variants(complex_variants):
    """Testing reading and parsing complex variants as CatalogScoreVariants"""
    pgs000021, pgs004239 = complex_variants

    # check no ValidationErrors are raised by the pydantic model, and:
    assert len(list(CatalogScoreVariant(**x) for x in pgs000021)) == 33
    assert len(list(CatalogScoreVariant(**x) for x in pgs004239)) == 14
