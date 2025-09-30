from __future__ import annotations

import zarr.codecs

# ./. (no call) needs to be represented in np.uint8 (np.nan requires float)
MISSING_GENOTYPE_SENTINEL_VALUE: int = 255
MISSING_GENOTYPE_TUPLE: tuple[int, int] = (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    MISSING_GENOTYPE_SENTINEL_VALUE,
)
# see bgen: 0 = AA / 1 = Aa / 2 = aa / 3 = missing
MISSING_GENOTYPE_FANCY_INDEX: int = 3

# the largest number of variants that can be stored in the zarr genotype array
# vibes-based choice: biggest PGS in 2025 are ~10,000,000 variants, so 10x that
ZARR_MAX_N_VARIANTS: int = 100_000_000
ZARR_PLOIDY: int = 2  # diploid

ZARR_VARIANT_CHUNK_SIZE: int = 1_000

# numerical arrays will compress quite well with lower compression levels
ZARR_COMPRESSOR = zarr.codecs.ZstdCodec(level=3)

BGEN_PHASED_N_COLS = 4
BGEN_UNPHASED_N_COLS = 3
