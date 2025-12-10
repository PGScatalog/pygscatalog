from __future__ import annotations

from typing import Literal

from numpy.dtypes import StringDType

# ./. (no call) needs to be represented in np.uint8 (np.nan requires float)
MISSING_GENOTYPE_SENTINEL_VALUE: Literal[255] = 255
MISSING_GENOTYPE_TUPLE = (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    MISSING_GENOTYPE_SENTINEL_VALUE,
)

ZARR_PLOIDY: Literal[2] = 2  # diploid

# variable-width string
# https://numpy.org/devdocs/user/basics.strings.html#variable-width-strings
NUMPY_STRING_DTYPE = StringDType()

# TODO: fix in next release  # noqa: TD002 TD003
VALID_CHROMOSOMES = {str(x) for x in range(1, 23)}.union({None})
