from __future__ import annotations

import logging
from typing import TYPE_CHECKING, TypedDict, cast

import dask.config
import numpy as np
import zarr
import zarr.storage
from dask import array as da

from pgscatalog.calc.lib.config import config
from pgscatalog.calc.lib.constants import (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    NO_ALLELE_SENTINEL_VALUE,
)

from ._impute import calculate_mean_dosage

if TYPE_CHECKING:
    import duckdb
    import polars as pl
    from numpy import typing as npt


logger = logging.getLogger(__name__)


def make_dosage_arrays(
    df: pl.DataFrame,
    zarr_group: str,
    genotypes_group: zarr.Group,
    pgs_store: zarr.storage.StoreLike,
    sampleset: str,
    n_samples: int,
    conn: duckdb.DuckDBPyConnection,
    n_workers: int = 1,
    n_minimum_impute: int = 50,
) -> tuple[zarr.Array, zarr.Array, zarr.Array]:
    """
    Calculate and adjust effect allele dosage, storing dosage and missingness masks in
    Zarr arrays (disk-backed).

    This function reads one shard's 3D genotype array to create a 2D dosage array
    where rows = samples and columns = dosage.

    Missing dosages are imputed from mean effect allele frequency.

    Dosage is adjusted for recessive and dominant effect types.

    Parameters
    ----------
    df : pl.DataFrame
        Score variants in the shard being scored, which must contain:
        - ``target_row_nr``: row indices of variants in the genotype array
        - ``effect_allele_idx``: allele index for effect allele dosage calculation
        - ``is_recessive``: boolean mask for recessive effect type
        - ``is_dominant``: boolean mask for dominant effect type
    zarr_group : str
        The shard being scored: its path relative to the store root
        (``sampleset/chromosome/filename``)
    genotypes_group: zarr.Group
        A read only group which contains genotype arrays
    pgs_store : zarr.storage.StoreLike
        A writable zarr store where dosage and missing arrays will be written
    sampleset : str
        Identifier of the sample set within the Zarr store
    n_samples : int
        Number of samples the shard holds, which the cache reports (it owns the
        hierarchy, so nothing here walks it)
    conn : duckdb.DuckDBPyConnection
        A connection to the score database, which holds the variant metadata
        (read-only is enough)
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
    ploidy_array : zarr.Array
        Uint8 array of shape ``(n_variants, n_samples)`` holding each call's observed
        ploidy, which is what ``allele_count`` sums

    Notes
    -----
    - Dosages are imputed as ``ploidy * mean(effect allele frequency)`` when
      sufficient non-missing calls are present
    - Both arrays are created with chunking along the variant dimension,
      using ``ZARR_VARIANT_CHUNK_SIZE`` and compressed with ``ZARR_COMPRESSOR``.
    """
    n_variants = get_n_score_variants(conn=conn)
    pgs_group = zarr.open_group(store=pgs_store, path=f"pgs/{sampleset}")

    # an array to record dosage of effect alleles
    dosage_array: zarr.Array = pgs_group.create_array(
        "dosage",
        shape=(n_variants, n_samples),
        dtype=np.float64,
        fill_value=np.nan,
        overwrite=True,
        chunks=(config.ZARR_VARIANT_CHUNK_SIZE, n_samples),
        compressors=config.ZARR_COMPRESSOR,
    )

    # an array to track which variants are missing
    missing_array: zarr.Array = pgs_group.create_array(
        "missing",
        shape=(n_variants, n_samples),
        dtype=np.bool_,
        fill_value=np.nan,
        overwrite=True,
        chunks=(config.ZARR_VARIANT_CHUNK_SIZE, n_samples),
        compressors=config.ZARR_COMPRESSOR,
    )

    # an array to record how many alleles each call actually holds, which is 2 for a
    # diploid call and 1 for a hemizygous one (ADR 0007)
    ploidy_array: zarr.Array = pgs_group.create_array(
        "ploidy",
        shape=(n_variants, n_samples),
        dtype=np.uint8,
        overwrite=True,
        chunks=(config.ZARR_VARIANT_CHUNK_SIZE, n_samples),
        compressors=config.ZARR_COMPRESSOR,
    )

    # calculate dosage for this shard's genotypes (sampleset/chromosome/filename) and
    # store it in the dosage array, which holds nothing else
    write_dosage_arrays(
        genotypes_group=genotypes_group,
        zarr_group=zarr_group,
        df=df,
        n_minimum_impute=n_minimum_impute,
        dosage_array=dosage_array,
        missing_array=missing_array,
        ploidy_array=ploidy_array,
        n_workers=n_workers,
    )

    # in these arrays the fill value (np.nan) has been replaced with computed values
    return missing_array, dosage_array, ploidy_array


def write_dosage_arrays(
    genotypes_group: zarr.Group,
    zarr_group: str,
    df: pl.DataFrame,
    n_minimum_impute: int,
    dosage_array: zarr.Array,
    missing_array: zarr.Array,
    ploidy_array: zarr.Array,
    n_workers: int = 1,
) -> None:
    """
    Compute one shard's genotype dosages and write them into Zarr arrays.

    For the shard being scored, this function:

    1. Loads genotype data from the genotypes (read-only) Zarr group.
    2. Calculates effect allele dosages and imputes missing values.
    3. Adjusts for dominant/recessive effect types.
    4. Writes the resulting arrays to the dosage and missingness arrays.

    Args:
        genotypes_group: Root Zarr group containing the shard.
        zarr_group: Path of the shard being scored, relative to the store root.
        df: DataFrame of the shard's score variants with variant metadata.
        n_minimum_impute: Minimum number of samples required to impute missing dosages.
        dosage_array: Target Zarr array to store computed dosages.
        missing_array: Target Zarr array to store missingness flags.
        ploidy_array: Target Zarr array to store each call's observed ploidy.

    Notes:
        The dosage_array and missing_array are initialised with np.nan. If any NaN
        values remain after the writes are finished an exception will be raised.

        The input arrays dosage_array and missing_array are modified in place.
    """
    logger.info(f"Calculating and storing dosage for {zarr_group=}")

    target_group = genotypes_group[zarr_group]

    if not isinstance(target_group, zarr.Group):
        raise TypeError

    target_zarr = target_group["genotypes"]

    target_idx = df["target_row_nr"].to_numpy().astype(np.int64)
    effect_idx = df["effect_allele_idx"].to_numpy()

    # extract score variants from cache
    sliced_gts: da.Array = da.from_zarr(
        target_zarr, chunks=(config.ZARR_VARIANT_CHUNK_SIZE, None, None)
    )[target_idx]

    # ploidy is derived from the genotype array, never stored in the cache (ADR 0007)
    ploidy: da.Array = da.sum(
        sliced_gts != NO_ALLELE_SENTINEL_VALUE, axis=2, dtype=np.uint8
    )

    # replace missing calls with np.nan
    missing_calls_masked: da.Array = sentinel_to_nan(sliced_gts)

    # calculate dosage of the effect allele
    dosage: da.Array = calculate_effect_allele_dosage(missing_calls_masked, effect_idx)

    # impute missing dosage: the call's own ploidy * mean effect allele frequency
    is_dosage_nan, filled = fill_missing_dosages(
        dosage_array=dosage, ploidy_array=ploidy, n_minimum_samples=n_minimum_impute
    )

    # adjust dosages for recessive / dominant effect types
    recessive_mask: npt.NDArray[np.bool_] = df["is_recessive"].to_numpy()
    dominant_mask: npt.NDArray[np.bool_] = df["is_dominant"].to_numpy()
    effect_type_adjusted: da.Array = adjust_dosage_for_effect(
        dosage_array=filled,
        recessive_mask=recessive_mask,
        dominant_mask=dominant_mask,
    )

    # the arrays hold this shard and nothing else, so each is written whole:
    # weight_mat_row_nr is the position of a row in this shard's ordering, so there
    # is no region to offset into
    dask_tasks = [
        effect_type_adjusted.to_zarr(url=dosage_array, compute=False),
        is_dosage_nan.to_zarr(url=missing_array, compute=False),
        ploidy.to_zarr(url=ploidy_array, compute=False),
    ]

    with dask.config.set(scheduler="threads", num_workers=n_workers):
        logger.info(
            f"Computing effect allele dosage and missing calls with {n_workers=}"
        )
        dask.compute(*dask_tasks)

        # the arrays were created full of NaN, so a write that didn't cover them
        # leaves NaNs behind - a wrong score, with nothing else to notice it
        def has_nan(zarr_array: zarr.Array) -> np.bool_:
            return cast(
                "np.bool_", da.any(da.isnan(da.from_zarr(zarr_array))).compute()
            )

        is_any_nan = has_nan(dosage_array) or has_nan(missing_array)

    logger.info("Finished calculating effect allele dosage")

    if is_any_nan:
        reason = no_called_alleles_reason(dosage_array, ploidy_array)
        logger.critical(
            "NaNs detected in dosage array or missing array after calculating "
            f"dosage: some variants have no called alleles because {reason}. "
            "Impossible to continue PGS calculation."
        )
        raise ValueError(
            "NaN values detected in dosage array or missing array: some variants "
            f"have no called alleles because {reason}"
        )
    logger.info("Dosage/missing arrays doesn't contain NaN values (good!)")


def no_called_alleles_reason(dosage_array: zarr.Array, ploidy_array: zarr.Array) -> str:
    """Why a variant had no called alleles to take an effect allele frequency from.

    There are two ways to get there and they need different messages: the second is
    biology, and reads as data corruption if you don't know that (ADR 0009).
    """
    nan_variants = da.any(da.isnan(da.from_zarr(dosage_array)), axis=1)
    nobody_can_carry = da.all(da.from_zarr(ploidy_array) == 0, axis=1)

    if bool(da.any(nan_variants & nobody_can_carry).compute()):
        return (
            "no sample can carry them (every call is nullploid, e.g. a chrY "
            "variant in a shard of women)"
        )
    return "every call is missing"


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

    Only the missing-call sentinel goes down this path. The no-allele sentinel
    must not: `calculate_effect_allele_dosage` NaNs a whole call if any of its
    slots is NaN, so every hemizygous call would be imputed instead of scored
    (ADR 0007). It needs no mask - no allele index equals it.

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
    ploidy_array: da.Array,
    n_minimum_samples: int = 50,
) -> tuple[da.Array, da.Array]:
    """
    Replace missing dosages with mean effect allele dosage.

    Missing dosage values (NaN) are replaced with the mean effect allele dosage
    calculated per (variant, sample), if at least `n_minimum_samples` non-missing
    values are available. Returns both a mask of originally missing values and the
    imputed dosage array.

    Parameters
    ----------
    dosage_array : dask.array.Array
        2D array of effect allele dosages with shape (n_variants, n_samples).
        Missing dosages are represented as NaN.
    ploidy_array : dask.array.Array
        2D array of the same shape holding each call's observed ploidy (ADR 0007)
    n_minimum_samples : int, default=50
        Minimum number of non-missing samples required to compute the mean dosage
        for imputation. An exception will be raised if too few samples are present
        *and* there is something to impute — a small shard with no missing calls
        needs no mean, so it is returned unchanged.

    Returns
    -------
    is_dosage_nan : dask.array.Array
        Boolean array of the same shape as `dosage_array` indicating which values
        were originally missing (NaN).
    filled : dask.array.Array
        Float64 array with missing dosages replaced by the mean for that call
    """
    is_dosage_nan = da.isnan(dosage_array)

    # ponytail: only shards below the floor pay for the missingness scan. A shard with
    # enough samples can't trip it, so don't walk the array to find that out.
    if dosage_array.shape[1] < n_minimum_samples and not bool(
        is_dosage_nan.any().compute()
    ):
        logger.info("No missing dosages to impute, skipping the sample floor check")
        return is_dosage_nan, dosage_array

    mean_effect_dosage: da.Array = calculate_mean_dosage(
        dosage_array=dosage_array,
        ploidy_array=ploidy_array,
        n_minimum_samples=n_minimum_samples,
    )

    filled: da.Array = da.where(is_dosage_nan, mean_effect_dosage, dosage_array)

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
        recessive_mask, chunks=(config.ZARR_VARIANT_CHUNK_SIZE,)
    )[:, np.newaxis]
    dominant_mask_2d = da.from_array(
        dominant_mask, chunks=(config.ZARR_VARIANT_CHUNK_SIZE,)
    )[:, np.newaxis]

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


def get_n_score_variants(conn: duckdb.DuckDBPyConnection) -> int:
    """
    Get the total number of score variants from the weight matrix.

    Parameters
    ----------
    conn : duckdb.DuckDBPyConnection
        A connection to the score database, which holds the `wide_score_variants`
        table (read-only is enough).

    Returns
    -------
    n_variants : int
        Number of variants in the `wide_score_variants` table.

    Raises
    ------
    ValueError
        If the query returns no results, which indicates bad database state.
    """
    query = conn.sql("SELECT COUNT(*) FROM wide_score_variants").fetchone()

    if query is None:
        logger.critical("Couldn't get number of variants from dosage table")
        raise ValueError("Query returned no results")
    return int(query[0])
