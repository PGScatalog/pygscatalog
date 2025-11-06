from __future__ import annotations

import logging
from typing import TYPE_CHECKING, TypedDict

import dask.config
import duckdb
import numpy as np
import zarr
import zarr.storage
from dask import array as da

from ..constants import (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    ZARR_COMPRESSOR,
    ZARR_VARIANT_CHUNK_SIZE,
)
from ._impute import calculate_mean_dosage

if TYPE_CHECKING:
    import polars as pl
    from numpy import typing as npt

    from ..types import Pathish


logger = logging.getLogger(__name__)


def store_dosage_from_chunks(
    df_groups: dict[str, pl.DataFrame],
    store: zarr.storage.StoreLike,
    sampleset: str,
    db_path: Pathish,
    n_workers: int = 1,
    n_minimum_impute: int = 50,
) -> tuple[zarr.Array, zarr.Array]:
    """
    Calculate and adjust effect allele dosage, storing dosage and missingness masks in
    Zarr arrays (disk-backed).

    This function basically queries multiple Zarr arrays in different groups to create
    a single 2D dosage array where rows = samples and columns = dosage.

    Missing dosages are imputed from mean effect allele frequency.

    Dosage is adjusted for recessive and dominant effect types.

    Parameters
    ----------
    df_groups : dict of str -> pl.DataFrame
        Mapping of Zarr genotype group names to dataframes describing score variants
        within that group. Each dataframe must contain:
        - ``target_row_nr``: row indices of variants in the genotype array
        - ``effect_allele_idx``: allele index for effect allele dosage calculation
        - ``is_recessive``: boolean mask for recessive effect type
        - ``is_dominant``: boolean mask for dominant effect type
    store : zarr.storage.StoreLike
        Zarr store where dosage and missing arrays will be written
    sampleset : str
        Identifier of the sample set within the Zarr store
    db_path : Pathish
        Path to the variant database which contains metadata
    n_workers : int, default=1
        Number of threads to use when writing dosage/missing arrays with Dask
    n_minimum_impute : int, default=50
        Minimum number of samples required for allele frequency-based imputation.
        If fewer are available an exception will be raised.

    Returns
    -------
    missing_array : zarr.Array
        Boolean array of shape ``(n_variants, n_samples)`` indicating missing/imputed
        dosage entries
    dosage_array : zarr.Array
        Float64 array of shape ``(n_variants, n_samples)`` storing effect allele dosages

    Notes
    -----
    - Dosages are imputed as ``2 * mean(effect allele frequency)`` when sufficient
      non-missing calls are present
    - Array writes are performed in blocks using the region argument to minimise I/O
    - Both arrays are created with chunking along the variant dimension,
      using ``ZARR_VARIANT_CHUNK_SIZE`` and compressed with ``ZARR_COMPRESSOR``.
    """
    n_samples = get_n_samples(store=store, sampleset=sampleset)
    n_variants = get_n_score_variants(db_path=db_path)
    pgs_group = zarr.open_group(store=store, path=f"pgs/{sampleset}")

    # an array to record dosage of effect alleles
    dosage_array: zarr.Array = pgs_group.create_array(
        "dosage",
        shape=(n_variants, n_samples),
        dtype=np.float64,
        fill_value=np.nan,
        overwrite=True,
        chunks=(ZARR_VARIANT_CHUNK_SIZE, n_samples),
        compressors=ZARR_COMPRESSOR,
    )

    # an array to track which variants are missing
    missing_array: zarr.Array = pgs_group.create_array(
        "missing",
        shape=(n_variants, n_samples),
        dtype=np.bool_,
        fill_value=np.nan,
        overwrite=True,
        chunks=(ZARR_VARIANT_CHUNK_SIZE, n_samples),
        compressors=ZARR_COMPRESSOR,
    )

    # prepare target gts
    gt_group = zarr.open_group(store=store, path="gts", mode="r")

    start = 0
    for zarr_group, df in df_groups.items():
        logger.info(f"Calculating and storing dosage for {zarr_group=}")
        target_group = gt_group[zarr_group]

        if not isinstance(target_group, zarr.Group):
            raise TypeError

        target_zarr = target_group["genotypes"]

        target_idx = df["target_row_nr"].to_numpy().astype(np.int64)
        effect_idx = df["effect_allele_idx"].to_numpy()

        # extract score variants from cache
        sliced_gts: da.Array = da.from_zarr(
            target_zarr, chunks=(ZARR_VARIANT_CHUNK_SIZE, None, None)
        )[target_idx]

        # replace missing calls with np.nan
        missing_calls_masked: da.Array = sentinel_to_nan(sliced_gts)

        # calculate dosage of the effect allele
        dosage: da.Array = calculate_effect_allele_dosage(
            missing_calls_masked, effect_idx
        )

        # impute missing dosage: 2 * mean effect allele frequency
        is_dosage_nan, filled = fill_missing_dosages(
            dosage_array=dosage, n_minimum_samples=n_minimum_impute
        )

        # adjust dosages for recessive / dominant effect types
        recessive_mask: npt.NDArray[np.bool_] = df["is_recessive"].to_numpy()
        dominant_mask: npt.NDArray[np.bool_] = df["is_dominant"].to_numpy()
        effect_type_adjusted: da.Array = adjust_dosage_for_effect(
            dosage_array=filled,
            recessive_mask=recessive_mask,
            dominant_mask=dominant_mask,
        )

        n_group_variants = target_idx.shape[0]
        end = start + n_group_variants

        # write specific regions to the zarr array to minimise memory usage
        with dask.config.set(scheduler="threads", num_workers=n_workers):
            logger.info(f"Writing dosage region {start=} {end=} to zarr array")
            effect_type_adjusted.to_zarr(
                url=dosage_array,
                component="dosage",
                # region is important to avoid writing the entire array
                region=(slice(start, end), slice(None)),
                compute=True,
            )
            logger.info("Finished writing dosage region")

            logger.info("Writing missing calls to array")
            is_dosage_nan.to_zarr(
                url=missing_array,
                component="missing",
                region=(slice(start, end), slice(None)),
                compute=True,
            )

        start = end
    logger.info("Finished calculating effect allele dosage")

    return missing_array, dosage_array


def calculate_effect_allele_dosage(
    genotype_array: da.Array, effect_idx: npt.NDArray[np.int64]
) -> da.Array:
    """
    Calculate effect allele dosages from a genotype call array.

    Each genotype call is compared against the effect allele index to determine
    how many copies of the effect allele are present per sample. Missing calls
    (NaN) propagate into the output dosage array.


    Parameters
    ----------
    genotype_array : dask.array.Array
        3D array of integer-encoded genotypes with shape
        ``(n_variants, n_samples, ploidy)``
        - Values correspond to allele indices
        - Missing calls are represented as ``NaN``
    effect_idx : numpy.ndarray of int64, shape (n_variants,)
        Index of the effect allele for each variant.

    Returns
    -------
    dosage : dask.array.Array
        2D float array of shape ``(n_variants, n_samples)`` giving the effect allele
        dosage per variant and sample. Values are:
        - ``0`` → no effect allele copies
        - ``1`` → heterozygous (one effect allele)
        - ``2`` → homozygous (two effect alleles)
        - ``NaN`` if any allele call for the sample/variant was missing.
    """
    # does each call match the effect allele index? e.g. effect allele == REF -> 0
    matches: da.Array = genotype_array == effect_idx[:, None, None]

    # but anything == np.nan will always be false, so find NaNs in the original array
    nan_mask: da.Array = da.isnan(genotype_array)

    # sum boolean matches along axis=2 to create a dosage array (n_variants, n_samples)
    dosage: da.Array = da.sum(matches, axis=2)

    # if any genotype in the ploidy slice was NaN, set the dosage to NaN
    masked_dosage: da.Array = da.where(da.any(nan_mask, axis=2), np.nan, dosage)

    return masked_dosage


def sentinel_to_nan(target_dask: da.Array) -> da.Array:
    """
    Replace sentinel genotype values with NaN and cast to float64.

    Genotype arrays may use a sentinel integer to represent missing values,
    which cannot be represented as NaN in integer arrays.

    Returns
    -------
    masked : dask.array.Array
        Float64 array with the same shape as `target_dask`, where sentinel values
        are replaced with NaN and all other values are cast to float64.
    """
    # integer arrays can't represent NaN
    floating_target = target_dask.astype(np.float64)
    sentinel = np.float32(MISSING_GENOTYPE_SENTINEL_VALUE)
    masked: da.Array = da.where(floating_target == sentinel, np.nan, floating_target)
    return masked


def fill_missing_dosages(
    dosage_array: da.Array,
    n_minimum_samples: int = 50,
) -> tuple[da.Array, da.Array]:
    """
    Replace missing dosages with mean effect allele dosage.

    Missing dosage values (NaN) are replaced with the mean effect allele dosage
    calculated per variant, if at least `n_minimum_samples` non-missing values
    are available. Returns both a mask of originally missing values and the
    imputed dosage array.

    Parameters
    ----------
    dosage_array : dask.array.Array
        2D array of effect allele dosages with shape (n_variants, n_samples).
        Missing dosages are represented as NaN.
    n_minimum_samples : int, default=50
        Minimum number of non-missing samples required to compute the mean dosage
        for imputation. An exception will be raised if too few samples are present.

    Returns
    -------
    is_dosage_nan : dask.array.Array
        Boolean array of the same shape as `dosage_array` indicating which values
        were originally missing (NaN).
    filled : dask.array.Array
        Float64 array with missing dosages replaced by the per-variant mean
    """
    is_dosage_nan = da.isnan(dosage_array)

    mean_effect_dosage: da.Array = calculate_mean_dosage(
        dosage_array=dosage_array,
        n_minimum_samples=n_minimum_samples,
    )

    # broadcast for per-row replacement
    mean_dosage = mean_effect_dosage[:, np.newaxis]

    filled: da.Array = da.where(is_dosage_nan, mean_dosage, dosage_array)

    return is_dosage_nan, filled


def adjust_dosage_for_effect(
    recessive_mask: npt.NDArray[np.bool_],
    dominant_mask: npt.NDArray[np.bool_],
    dosage_array: da.Array,
) -> da.Array:
    """
    Adjust effect allele dosages according to variant effect types.

    Dosages are modified based on the variant's effect type:
      - Additive variants are unchanged
      - Dominant variants are capped at 1 (any dosage > 1 becomes 1)
      - Recessive variants are reduced by 1 and clamped at 0
        (i.e., dosage = max(dosage - 1, 0))

    Parameters
    ----------
    recessive_mask : numpy.ndarray of bool, shape (n_variants,)
        Boolean mask indicating which variants have a recessive effect type.
    dominant_mask : numpy.ndarray of bool, shape (n_variants,)
        Boolean mask indicating which variants have a dominant effect type.
    dosage_array : dask.array.Array
        2D array of effect allele dosages with shape (n_variants, n_samples).

    Returns
    -------
    adjusted_dosage : dask.array.Array
        Dosage array with the same shape as `dosage_array`, with adjustments
        applied for recessive and dominant variants.

    Raises
    ------
    ValueError
        If any variant is marked as both dominant and recessive.
    """
    if np.sum(recessive_mask) == 0 and np.sum(dominant_mask) == 0:
        logger.info("No variants need to be adjusted for effect type, skipping")
        return dosage_array
    logger.info("Adjusting effect allele dosage for dominant/recessive effect types")

    if np.any(recessive_mask & dominant_mask):
        raise ValueError("A variant cannot be both dominant and recessive")

    # expand 1D masks to 2D
    recessive_mask_2d = da.from_array(
        recessive_mask, chunks=(ZARR_VARIANT_CHUNK_SIZE,)
    )[:, np.newaxis]
    dominant_mask_2d = da.from_array(dominant_mask, chunks=(ZARR_VARIANT_CHUNK_SIZE,))[
        :, np.newaxis
    ]

    # apply dominant adjustment (cap at 1)
    dominant_adjusted = da.where(
        dominant_mask_2d, da.minimum(dosage_array, 1), dosage_array
    )

    # apply recessive adjustment (subtract 1 and clamp to 0)
    recessive_adjusted: da.Array = da.where(
        recessive_mask_2d, da.maximum(dominant_adjusted - 1, 0), dominant_adjusted
    )

    return recessive_adjusted


class ScoreStats(TypedDict):
    """Container for score statistics used to build a Polars DataFrame."""

    # per-sampleset data
    sampleset: list[str]
    sample_id: list[list[str]]

    accession: list[str]
    # per-sample weighted sum for each accession
    score: list[npt.NDArray[np.float64]]
    # per-sample dosage sums for each accession
    dosage_sum: list[npt.NDArray[np.float64]]
    # per-sample: non-missing allele counts for each accession
    allele_count: list[npt.NDArray[np.int64]]
    # per-accession: number of matched variants
    n_matched: list[int]


def get_n_score_variants(db_path: Pathish) -> int:
    """
    Get the total number of score variants from the weight matrix.

    Parameters
    ----------
    db_path : Pathish
        Path to the DuckDB database containing the `wide_score_variants` table.

    Returns
    -------
    n_variants : int
        Number of variants in the `wide_score_variants` table.

    Raises
    ------
    ValueError
        If the query returns no results, which indicates bad database state.
    """
    with duckdb.connect(str(db_path), read_only=True) as conn:
        query = conn.sql("SELECT COUNT(*) FROM wide_score_variants").fetchone()

    if query is None:
        logger.critical("Couldn't get number of variants from dosage table")
        raise ValueError("Query returned no results")
    return int(query[0])


def get_n_samples(store: zarr.storage.StoreLike, sampleset: str) -> int:
    """
    Get the number of samples in a given sample set stored in Zarr.

    Parameters
    ----------
    store : zarr.storage.StoreLike
        Zarr store containing genotype data.
    sampleset : str
        Name of the sample set (group) in the store.

    Returns
    -------
    n_samples : int
        Number of samples in the sample set.

    Raises
    ------
    ValueError
        If the `samples` attribute in the Zarr group is not a list.
    """
    group = zarr.open_group(store=store, path=f"gts/{sampleset}", mode="r")
    samples = group.attrs["samples"]
    if not isinstance(samples, list):
        raise ValueError("Expected samples attribute to be a list")

    return len(samples)
