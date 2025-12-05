from __future__ import annotations

import logging

import numpy as np
from dask import array as da

from pgscatalog.calc.lib.constants import ZARR_PLOIDY

logger = logging.getLogger(__name__)


def calculate_mean_dosage(
    dosage_array: da.Array, 
    n_minimum_samples: int = 50
    ) -> da.Array:
    """
    Impute effect allele dosage using the allelic frequency of the effect allele

    Parameters
    ----------
    dosage_array : da.Array
        Shape (n_variants, n_samples): np.nan represents missing calls
    n_workers : int
        Number of threads for dask to use
    n_minimum_samples : int
        Minimum number of samples required to impute using allelic frequency.

    Raises a ValueError if the genotype array has fewer than n_minimum_samples samples.

    Returns
    -------
    npt.NDArray[np.float64]
        Shape (variants, ), the imputed effect allele dosage
    """
    if not np.issubdtype(dosage_array.dtype, np.floating):
        raise TypeError("The genotype array must contain floats")

    if (n_samples := dosage_array.shape[1]) < n_minimum_samples:
        raise ValueError(f"{n_samples=} is less than {n_minimum_samples=}")

    logger.info("Calculating mean dosage from allelic frequency")

    # numerator: sum of effect alleles per variant, treating NaN as 0
    n_effect_alleles = da.nansum(dosage_array, axis=1)  # shape (n_variants,)

    # per-variant number of **called genotypes**
    n_called_genotypes = da.sum(~da.isnan(dosage_array), axis=1)  # shape (n_variants,)

    # denominator: number of called alleles per variant (called genotypes * ploidy)
    # replace this with a ploidy mask when X / Y / MT support added
    n_called_alleles = n_called_genotypes * ZARR_PLOIDY

    # where no called alleles, set freq to NaN
    # scores get checked for NaN values before writing, so that's OK
    effect_allele_freq = da.where(
        n_called_alleles > 0,
        n_effect_alleles / n_called_alleles,
        da.array(np.nan, dtype=float),
    )

    mean_effect_allele_dosage = effect_allele_freq * ZARR_PLOIDY

    return mean_effect_allele_dosage