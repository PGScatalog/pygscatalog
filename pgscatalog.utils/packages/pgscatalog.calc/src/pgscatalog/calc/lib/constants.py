from __future__ import annotations

from functools import cache
from typing import Literal, overload

from numpy.dtypes import StringDType

# ./. (no call) needs to be represented in np.uint8 (np.nan requires float)
MISSING_GENOTYPE_SENTINEL_VALUE: Literal[255] = 255
MISSING_GENOTYPE_TUPLE = (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    MISSING_GENOTYPE_SENTINEL_VALUE,
)

# an allele slot the call doesn't have: a hemizygous call is [allele, 254]. Distinct
# from a missing call, which is a slot the sample has but the caller couldn't fill.
# Ploidy is sum(gt != 254) and is derived, never stored (ADR 0007)
NO_ALLELE_SENTINEL_VALUE: Literal[254] = 254

# the genotype axis is as wide as the widest call, not the ploidy of every call
ZARR_MAX_PLOIDY: Literal[2] = 2

# variable-width string
# https://numpy.org/devdocs/user/basics.strings.html#variable-width-strings
NUMPY_STRING_DTYPE = StringDType()

# the LOAD gate, now the whole canonical vocabulary: X and Y read haploid calls and
# MT collapses diploid-coded ones (ADR 0007). A mixed-sex chrY VCF still gets no
# further than the VCF reader, which can't tell a `.` from a missing Y
VALID_CHROMOSOMES = {str(x) for x in range(1, 23)}.union({"X", "Y", "MT", None})

# characters that can't appear in a plain base sequence, and so mark out the VCF
# alleles that name no sequence at all: the spanning deletion (*), symbolic alleles
# (<NON_REF>, <DEL>, ...), breakends (]13:123]T, T[2:456[, .TAGT) and the missing
# allele (.). This is what separates a legal allele no score can match from a corrupt
# one: XYZ contains none of these, so it stays a ValueError.
_NON_SEQUENCE_CHARACTERS: frozenset[str] = frozenset("*<>[].")


@cache
def is_unmatchable_allele(allele: str) -> bool:
    """Is this a legal VCF allele that names no sequence to match against?

    Such alleles are kept in place and skipped at match time, never filtered away:
    a genotype call names its alleles by index, so dropping one renames every allele
    after it (docs/adr/0010).

    Shared by the cache layer, which stores such an allele as written instead of
    refusing the record, and by the matcher, which skips it. Both layers import it
    from here so they agree on what is unmatchable by construction: the cache is the
    portable interface below scoring, so it cannot import the rule from the matcher.

    >>> is_unmatchable_allele("*")  # spanning deletion
    True
    >>> is_unmatchable_allele("<NON_REF>")  # symbolic allele
    True
    >>> is_unmatchable_allele("A"), is_unmatchable_allele("GATC")
    (False, False)

    Corrupt allele strings aren't unmatchable, they're wrong, and are left to raise:

    >>> is_unmatchable_allele("XYZ")
    False
    """
    return not _NON_SEQUENCE_CHARACTERS.isdisjoint(allele)


# dialects that survive stripping a chr prefix and zero padding. XY collapses into X
# because pseudoautosomal variants are ordinary X variants (ADR 0008), and 23-26 are
# plink's numeric codes
_CHROM_ALIASES = {"XY": "X", "M": "MT", "23": "X", "24": "Y", "25": "X", "26": "MT"}


@overload
def canonicalise_chrom(chrom: str | int) -> str: ...


@overload
def canonicalise_chrom(chrom: None) -> None: ...


@cache
def canonicalise_chrom(chrom: str | int | None) -> str | None:
    """Translate an input dialect's chromosome name into the internal vocabulary
    (`1`-`22`, `X`, `Y`, `MT`), so no stage below LOAD compares names across dialects.

    Called once at each boundary - every GenomeFileHandler on the way in, and the
    scoring file loader on the way in (ADR 0008).

    Parameters
    ----------
    chrom : str | int | None
        A chromosome name in any dialect, or None for unsplit whole-genome input,
        which stays None (`leaf_group_path` turns it into the `ALL` sentinel).
    """
    if chrom is None:
        return None

    s = str(chrom).removeprefix("chr").lstrip("0")
    return _CHROM_ALIASES.get(s, s)
