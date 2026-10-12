import pytest

from pgscatalog.calc.lib.constants import VALID_CHROMOSOMES, canonicalise_chrom


@pytest.mark.parametrize(
    "chrom, expected",
    [
        # scoring files
        ("X", "X"),
        ("XY", "X"),  # pseudoautosomal variants are ordinary X variants
        ("MT", "MT"),
        # VCF
        ("chr1", "1"),
        ("chrX", "X"),
        ("chrM", "MT"),
        ("M", "MT"),
        # bgen
        ("01", "1"),
        ("0X", "X"),
        ("chr01", "1"),
        # plink numeric codes
        ("23", "X"),
        ("24", "Y"),
        ("25", "X"),
        ("26", "MT"),
        # already canonical, and ints
        ("22", "22"),
        (7, "7"),
        # unsplit whole-genome input, which leaf_group_path turns into ALL
        (None, None),
    ],
)
def test_canonicalise_chrom(chrom, expected):
    assert canonicalise_chrom(chrom) == expected


def test_load_gate_is_open_for_the_whole_canonical_vocabulary():
    """X and Y read haploid calls and MT collapses diploid-coded ones (ADR 0007),
    so the gate and the vocabulary are finally the same set. Anything outside it -
    an unplaced contig, say - is still not a chromosome this pipeline scores"""
    assert {str(x) for x in range(1, 23)} | {"X", "Y", "MT", None} == VALID_CHROMOSOMES
    assert canonicalise_chrom("chrUn_KI270302v1") not in VALID_CHROMOSOMES
