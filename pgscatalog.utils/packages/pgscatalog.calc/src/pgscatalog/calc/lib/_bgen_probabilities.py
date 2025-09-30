import logging

import dask.array as da
import numpy as np
import numpy.typing as npt

from .constants import (
    BGEN_PHASED_N_COLS,
    BGEN_UNPHASED_N_COLS,
    MISSING_GENOTYPE_FANCY_INDEX,
    MISSING_GENOTYPE_SENTINEL_VALUE,
    ZARR_VARIANT_CHUNK_SIZE,
)

logger = logging.getLogger(__name__)


def unphased_probabilities_to_hard_calls(
    probabilities: list[npt.NDArray[np.floating]], n_samples: int
) -> da.Array:
    """Convert unphased probabilities to VCF style genotypes"""
    if not all(x.shape[1] == BGEN_UNPHASED_N_COLS for x in probabilities):
        raise ValueError("Unphased probabilities must have 3 columns")
    else:
        logger.info("Unphased probabilities have good shapes")

    # axis 0 (x): variants
    # axis 1 (y): samples
    # axis 2 (z): ploidy
    darr = da.from_array(
        np.stack(probabilities, axis=0),
        chunks=(ZARR_VARIANT_CHUNK_SIZE, n_samples, BGEN_UNPHASED_N_COLS),
    )

    # which variants have missing calls?
    missing_mask = da.any(darr == MISSING_GENOTYPE_SENTINEL_VALUE, axis=2)
    # which index is the most probable?
    allele_indices = darr.argmax(axis=2).astype(np.uint8)
    # reset the sentinel value to the fancy indexing missing value
    masked_indices = da.where(
        missing_mask, MISSING_GENOTYPE_FANCY_INDEX, allele_indices
    )

    # fancy indexing with a lookup table
    # 0: homozygous first allele
    # 1: heterozygous
    # 2: homozygous second allele
    # 3: missing
    # (lookup table is converting maximum probabilities to hard calls)
    lookup_table = np.array(
        [
            [0, 0],
            [0, 1],
            [1, 1],
            [MISSING_GENOTYPE_SENTINEL_VALUE, MISSING_GENOTYPE_SENTINEL_VALUE],
        ],
        dtype=np.uint8,
    )

    def map_lookup(block):
        return lookup_table[block]

    # dask doesn't support native fancy indexing, so map the numpy approach
    hard_calls = masked_indices.map_blocks(
        map_lookup,
        dtype=np.uint8,
        new_axis=2,
    )

    return hard_calls


def phased_probabilities_to_hard_calls(
    probabilities: list[npt.NDArray[np.floating]], n_samples: int
) -> da.Array:
    """Convert phased probabilities to genotypes"""
    if not all(x.shape[1] == BGEN_PHASED_N_COLS for x in probabilities):
        raise ValueError("Haplotype pair probabilities required (4 columns).")
    else:
        logger.info("Haplotype pair probabilities have good shapes")
        logger.info("Hard calling phased probabilities with dask")

    # axis 0 (x): variants
    # axis 1 (y): samples
    # axis 2 (z): ploidy
    darr = da.from_array(
        np.stack(probabilities, axis=0),
        chunks=(ZARR_VARIANT_CHUNK_SIZE, n_samples, BGEN_PHASED_N_COLS),
    )

    hap1, hap2 = darr[:, :, :2], darr[:, :, 2:]

    # convert to hard call: which index has the highest probability?
    hap1_alleles = da.argmax(hap1, axis=2).astype(np.uint8)
    hap2_alleles = da.argmax(hap2, axis=2).astype(np.uint8)

    # argmax will have broken the sentinel (missing) values
    # (because it returns the index of the sentinel value, which is large)
    hap1_mask = da.any(hap1 == MISSING_GENOTYPE_SENTINEL_VALUE, axis=2)
    hap2_mask = da.any(hap2 == MISSING_GENOTYPE_SENTINEL_VALUE, axis=2)

    # reset sentinel values in haplotypes
    masked_hap1 = da.where(hap1_mask, MISSING_GENOTYPE_SENTINEL_VALUE, hap1_alleles)
    masked_hap2 = da.where(hap2_mask, MISSING_GENOTYPE_SENTINEL_VALUE, hap2_alleles)

    return da.stack([masked_hap1, masked_hap2], axis=2)
