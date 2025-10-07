from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import duckdb
import numpy as np
import polars as pl
from dask import array as da
from numpy import typing as npt

from .constants import ZARR_VARIANT_CHUNK_SIZE

if TYPE_CHECKING:
    import zarr

    from ._dosage import ScoreStats
    from .types import Pathish


logger = logging.getLogger(__name__)


def create_score_table(db_path: Pathish) -> None:
    """Create the score table, which contains one row per sample per score

    Parameters
    ----------
    db_path : Pathish
        Path to the score database
    """
    logger.info("Creating or replacing score table")
    with duckdb.connect(str(db_path)) as conn:
        conn.sql("""
        CREATE OR REPLACE TABLE score_table (
            sampleset TEXT NOT NULL,
            accession TEXT NOT NULL,
            n_matched UINTEGER NOT NULL,
            sample_id TEXT NOT NULL,
            allele_count UINTEGER NOT NULL,
            score DOUBLE NOT NULL,
            dosage_sum DOUBLE NOT NULL,
            score_avg GENERATED ALWAYS AS (dosage_sum / allele_count) VIRTUAL,
            PRIMARY KEY (accession, sampleset, sample_id)
        );
        """)


def insert_scores(db_path: Pathish, _score_df: pl.DataFrame) -> None:
    # _ used because duckdb scans the polars dataframe variable but linters complain
    with duckdb.connect(str(db_path)) as conn:
        conn.sql("""
                 INSERT INTO score_table
                 SELECT sampleset,
                        accession,
                        n_matched,
                        sample_id,
                        allele_count,
                        score,
                        dosage_sum
                 FROM _score_df
                 """)


def write_scores(
    db_path: Pathish, out_dir: Pathish, max_memory_gb: str, threads: int, scale: int = 6
) -> None:
    """Write the score table to a directory in hive-partitioned format

    Scores are written out with a fixed width to preserve privacy (rounding limits
    identifiability). All calculations are done using internal cache data which is
    stored using full precision.

    Parameters
    ----------
    db_path : Pathish
        Path to the score database
    out_dir : Pathish
        Path to the output directory. Files in the directory will be overwritten.
    max_memory_gb : str
        Maximum RAM that duckdb can use before intermediate results spilling to disk
    threads : int
        Maximum number of threads duckdb can use
    scale : int
        Number of digits after the decimal point for score columns
    """
    logger.info(f"Copying score table to {out_dir}")
    logger.info(f"Score columns will have {scale} digits after the decimal point")
    with duckdb.connect(
        str(db_path),
        config={"max_memory": max_memory_gb, "threads": str(threads)},
    ) as conn:
        conn.sql(f"""
        COPY (
            SELECT
                sampleset,
                accession,
                sample_id,
                n_matched,
                allele_count,
                -- export fixed-width decimals
                CAST(dosage_sum AS DECIMAL(18, {scale})) AS dosage_sum,
                CAST(score AS DECIMAL(18, {scale})) AS score,
                CAST(score_avg AS DECIMAL(18, {scale})) AS score_avg,
            FROM score_table
            ORDER BY sampleset, accession, sample_id
        )
        TO '{str(out_dir)}'
        (
            FORMAT CSV,
            HEADER,
            DELIMITER ',',
            COMPRESSION GZIP,
            PARTITION_BY (sampleset, accession),
            WRITE_PARTITION_COLUMNS true,
            OVERWRITE true
        );
        """)
    logger.info("Finished copying score table")


def calculate_scores(dosage_array: zarr.Array, effect_weights: zarr.Array) -> da.Array:
    """
    Calculate a polygenic score

    Parameters
    ----------
    dosage_array : zarr.Array
        2D array of shape (variants, samples) containing effect allele dosage values
    effect_weights : zarr.Array
        1D array of shape (variants, ) containing effect sizes
    n_workers : int
        Number of threads to use for dask computation

    Returns
    -------
    np.ndarray
        1D NumPy array of scores for each sample
    """
    dosage = da.from_zarr(dosage_array, chunks=(ZARR_VARIANT_CHUNK_SIZE, None))
    weights = da.from_zarr(effect_weights, chunks=(ZARR_VARIANT_CHUNK_SIZE, None))
    if not dosage.shape[0] == weights.shape[0]:
        raise ValueError("dosage and weights must have same shape in axis 0")

    logger.info("Calculating weighted sum")
    logger.debug(f"{dosage.shape=}, {dosage.dtype=}")
    logger.debug(f"{weights.shape=}, {weights.dtype=}")

    scores: da.Array = da.dot(weights.T, dosage).T
    return scores


def calculate_score_statistics(
    db_path: Pathish,
    accessions: list[str],
    dosage_array: zarr.Array,
    is_missing_array: zarr.Array,
    sampleset: str,
    score: npt.NDArray[np.float64],
    sample_ids: list[str],
) -> pl.DataFrame:
    """
    Compute dosage statistics for a set of accessions.

    For each accession, the function:

      - Retrieves the number of matched variants from the allele match table
      - Selects variants with non-zero effect weights for that accession
      - Computes the sum of effect allele dosages across variants (per sample)
      - Computes the number of non-missing calls across variants (per sample)

    Parameters
    ----------
    db_path : Pathish
        Path to the DuckDB database file.
    sampleset : str
        Identifier for the sample set in the allele match table.
    accessions : list[str]
        List of accession identifiers to compute statistics for.
    dosage_array : zarr.Array
        2D dosage matrix with shape (variants, samples).

    Returns
    -------
    Polars DataFrame
        DataFrame with the following columns:
          - accession: score identifier
          - dosage_sum: sum of effect allele dosages
          - allele_count: number of variant sites with non-missing calls
          - n_matched: number of matched variants (constant per accession)
          - sample_id: sample identifier
          - score: weighted sum of effect allele dosages
          - sampleset: sampleset identifier

    This structure is important to cast to a polars dataframe.

    Notes
    -----
    To grab variants in an accession any rows with an effect weight of 0 are filtered
    from the accession column. If an accession contains variants with effect weights
    of zero (which would be very weird), these would get filtered out before calculating
    these statistics.
    """
    logger.info("Calculating dosage statistics")

    score_stats: ScoreStats = {
        "accession": [],
        "dosage_sum": [],
        "allele_count": [],
        "n_matched": [],
        "sample_id": [],
        "score": [],
        "sampleset": [],
    }

    with duckdb.connect(str(db_path), read_only=True) as conn:
        for accession in accessions:
            query = conn.execute(
                """
            SELECT COUNT(*)
            FROM allele_match_table
            WHERE match_result.match_summary='matched' AND
                sampleset = $sampleset AND
                accession = $accession;
            """,
                {"accession": accession, "sampleset": sampleset},
            ).fetchone()

            if query is None:
                raise ValueError(f"{accession=} can't count matches")
            n_matched = query[0]

            # slice the weight matrix to include variants with non-zero effect weights
            accession_idx = conn.sql(f"""
            SELECT weight_mat_row_nr
            FROM wide_score_variants
            WHERE "{accession}" != 0
            """).fetchnumpy()["weight_mat_row_nr"]

            accession_dosage = calculate_effect_allele_dosage(
                dosage_array=dosage_array, accession_idx=accession_idx
            )
            allele_count = calculate_nonmissing_allele_count(
                is_missing_array=is_missing_array,
                dosage_array=dosage_array,
                accession_idx=accession_idx,
            )

            # per-accession statistics
            score_stats["accession"].append(accession)
            score_stats["dosage_sum"].append(accession_dosage)
            score_stats["allele_count"].append(allele_count)
            score_stats["n_matched"].append(n_matched)

    # per-sampleset statistics
    score_stats["score"] = [score for score in score.T]  # noqa: C416
    score_stats["sampleset"] = [sampleset] * len(accessions)
    score_stats["sample_id"] = [sample_ids] * len(accessions)

    return pl.DataFrame(score_stats).explode(
        ["sample_id", "dosage_sum", "allele_count", "score"]
    )


def calculate_effect_allele_dosage(
    dosage_array: zarr.Array, accession_idx: npt.NDArray[np.int64]
) -> npt.NDArray[np.float64]:
    dosage_dask: da.Array = da.from_zarr(
        dosage_array, chunks=(ZARR_VARIANT_CHUNK_SIZE, None)
    )
    # returns a per-sample sum of dosages for all variants in this accession
    per_sample_dosage: npt.NDArray[np.float64] = (
        dosage_dask[accession_idx,].sum(axis=0).compute()
    )
    return per_sample_dosage


def calculate_nonmissing_allele_count(
    is_missing_array: zarr.Array,
    dosage_array: zarr.Array,
    accession_idx: npt.NDArray[np.int64],
) -> npt.NDArray[np.int64]:
    dosage_dask: da.Array = da.from_zarr(
        dosage_array, chunks=(ZARR_VARIANT_CHUNK_SIZE, None)
    )
    missing_dask: da.Array = da.from_zarr(
        is_missing_array, chunks=(ZARR_VARIANT_CHUNK_SIZE, None)
    )

    # missing variants will have been imputed, so replace these positions with np.nan
    missing_variants = da.where(~missing_dask, dosage_dask, np.nan)

    # return a per-sample sum of non-missing allele count in this accession
    per_sample_allele_count: npt.NDArray[np.int64] = da.sum(
        ~da.isnan(missing_variants[accession_idx,]), axis=0
    ).compute()
    return per_sample_allele_count
