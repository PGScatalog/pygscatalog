import pytest

from pgscatalog.calc.lib._bgen import add_chrom_padding, add_chrom_prefix


@pytest.mark.parametrize(
    "input_chrom, expected",
    [
        ("1", "chr1"),
        ("22", "chr22"),
        ("X", "chrX"),
    ],
)
def test_chrom_prefix(input_chrom, expected):
    assert add_chrom_prefix(input_chrom) == expected


@pytest.mark.parametrize(
    "input_chrom, expected",
    [
        ("1", "01"),
        ("11", "11"),
        ("X", "X"),
    ],
)
def test_chrom_padding(input_chrom, expected):
    assert add_chrom_padding(input_chrom) == expected
