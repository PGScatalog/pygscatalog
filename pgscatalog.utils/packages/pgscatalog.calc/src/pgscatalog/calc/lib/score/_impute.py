from __future__ import annotations

import logging

import numpy as np
from dask import array as da

from pgscatalog.calc.lib.constants import ZARR_PLOIDY

logger = logging.getLogger(__name__)


def calculate_mean_dosage(
    dosage_array: da.Array, n_minimum_samples: int = 50
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
    if not np.issubdtype(dosage_array.dtype, np.float64):
        raise TypeError("The genotype array must contain floats")

    if (n_samples := dosage_array.shape[1]) < n_minimum_samples:
        raise ValueError(f"{n_samples=} is less than {n_minimum_samples=}")

    logger.info("Calculating mean dosage from allelic frequency")

    # count effect alleles per variant
    n_effect_alleles = da.nansum(dosage_array, axis=1)

    # get number of alleles called in total
    n_alleles = dosage_array.size * ZARR_PLOIDY

    # allelic frequency = count of effect allele / number of called alleles
    effect_allele_freq = n_effect_alleles / n_alleles

    # mean dosage = effect allele frequency * 2
    mean_effect_dosage: da.Array = effect_allele_freq * ZARR_PLOIDY

    return mean_effect_dosage
