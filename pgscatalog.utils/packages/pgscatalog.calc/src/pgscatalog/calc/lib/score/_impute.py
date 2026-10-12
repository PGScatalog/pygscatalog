from __future__ import annotations

import logging

import numpy as np
from dask import array as da

logger = logging.getLogger(__name__)


def calculate_mean_dosage(
    dosage_array: da.Array, ploidy_array: da.Array, n_minimum_samples: int = 50
) -> da.Array:
    """
    Impute effect allele dosage using the allelic frequency of the effect allele

    Parameters
    ----------
    dosage_array : da.Array
        Shape (n_variants, n_samples): np.nan represents missing calls
    ploidy_array : da.Array
        Shape (n_variants, n_samples): how many alleles each call holds (ADR 0007)
    n_minimum_samples : int
        Minimum number of samples required to impute using allelic frequency.

    Raises a ValueError if the genotype array has fewer than n_minimum_samples samples.

    Returns
    -------
    npt.NDArray[np.float64]
        Shape (n_variants, n_samples), the imputed effect allele dosage
    """
    if not np.issubdtype(dosage_array.dtype, np.floating):
        raise TypeError("The genotype array must contain floats")

    if (n_samples := dosage_array.shape[1]) < n_minimum_samples:
        raise ValueError(f"{n_samples=} is less than {n_minimum_samples=}")

    logger.info("Calculating mean dosage from allelic frequency")

    # numerator: sum of effect alleles per variant, treating NaN as 0
    n_effect_alleles = da.nansum(dosage_array, axis=1)  # shape (n_variants,)

    # denominator: alleles actually called per variant, which is each call's own
    # ploidy rather than a constant. The frequency is pooled across ploidy - total
    # effect alleles over total called alleles - so hemizygous samples don't halve
    # it, and nullploid samples drop out for free (ADR 0009)
    n_called_alleles = da.sum(
        da.where(da.isnan(dosage_array), 0, ploidy_array), axis=1
    )  # shape (n_variants,)

    # where no called alleles, the frequency comes out NaN (dividing by NaN rather
    # than by 0 to keep numpy quiet about it)
    # scores get checked for NaN values before writing, so that's OK
    effect_allele_freq = n_effect_alleles / da.where(
        n_called_alleles > 0, n_called_alleles, np.nan
    )

    # ...but the imputed dosage is per sample: a man's missing X call is one expected
    # copy of the effect allele, a woman's is two
    mean_effect_allele_dosage: da.Array = ploidy_array * effect_allele_freq[:, None]
    return mean_effect_allele_dosage
