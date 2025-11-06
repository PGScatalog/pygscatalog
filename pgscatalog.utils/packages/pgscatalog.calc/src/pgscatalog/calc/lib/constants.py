from __future__ import annotations

from typing import Literal

import zarr.codecs

# ./. (no call) needs to be represented in np.uint8 (np.nan requires float)
MISSING_GENOTYPE_SENTINEL_VALUE: Literal[255] = 255
MISSING_GENOTYPE_TUPLE = (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    MISSING_GENOTYPE_SENTINEL_VALUE,
)

ZARR_PLOIDY: int = 2  # diploid

ZARR_VARIANT_CHUNK_SIZE: int = 5_000

# numerical arrays will compress quite well with lower compression levels
ZARR_COMPRESSOR = zarr.codecs.ZstdCodec(level=5)

# variant metadata is cached in fixed with string arrays
NUMPY_STRING_DTYPE = "<U255"

# TODO: fix in next release
VALID_CHROMOSOMES = {str(x) for x in range(1, 23)}.union({None})
