from __future__ import annotations

from functools import cache
from typing import Literal

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

# dialects that survive stripping a chr prefix and zero padding. XY collapses into X
# because pseudoautosomal variants are ordinary X variants (ADR 0008), and 23-26 are
# plink's numeric codes
_CHROM_ALIASES = {"XY": "X", "M": "MT", "23": "X", "24": "Y", "25": "X", "26": "MT"}


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
