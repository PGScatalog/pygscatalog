"""SCORE: calculate partial polygenic scores for one shard.

Staged with one shard, because `weight_mat_row_nr` indexes this task's dosage array
and staging two at once would silently misalign it. The dosage and weight arrays
it builds are scratch, not output: only `partial_scores.parquet` is published.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, cast

import dask.config
import zarr
import zarr.storage

from pgscatalog.calc.lib.cache.genotypecache import GenotypeCache

from . import stage
from ._dosage import make_dosage_arrays
from ._pgs import (
    calculate_score_statistics,
    calculate_scores,
    create_score_table,
    insert_scores,
)
from ._weight_matrix import store_group_weight_arrays
from .scoredb import ScoreDatabase

if TYPE_CHECKING:
    from collections.abc import Sequence

    import numpy as np
    import numpy.typing as npt
    import polars as pl

    from pgscatalog.calc.lib.types import Pathish, PathishList


logger = logging.getLogger(__name__)


class ScoreStage:
    """One SCORE task: one shard's genotypes, and the arrays derived from them.

    Parameters
    ----------
    out_dir : Pathish
        The task's output directory. `partial_scores.parquet` is published at its
        root; the cache, the database and the dosage arrays live in its scratch.
    max_memory_gb : float
        Maximum RAM duckdb may use before intermediate results spill to disk.
    threads : int
        Maximum number of threads duckdb and dask may use.
    minimum_samples_for_impute : int
        Minimum number of samples required to impute a missing dosage from mean
        effect allele frequency. Only SCORE calculates dosage, so only SCORE has
        anything to say about the floor.
    """

    def __init__(
        self,
        out_dir: Pathish,
        max_memory_gb: float,
        threads: int,
        minimum_samples_for_impute: int = 50,
    ):
        self._scratch = stage.scratch_dir(out_dir)

        # the cache owns the zarr hierarchy, so this stage asks it what it holds
        # instead of deriving paths. SCORE is staged with a cache, so it reads
        # the sampleset name and the shard from it
        self._cache = GenotypeCache(path=self._scratch / "genotypes.zarr")

        self._n_minimum_samples_for_impute = minimum_samples_for_impute

        # every query this stage runs goes through the database, so the resource
        # limits the CLI was given reach all of them
        self._db = ScoreDatabase(
            path=self._scratch / "pgs.db",
            max_memory_gb=max_memory_gb,
            threads=threads,
        )
        self._threads = threads

    def load_scores(self, scorefile_paths: Pathish | PathishList) -> None:
        """Copy scores processed with pgscatalog-format into a duckDB database table
        'scorevariants'
        """
        stage.load_scores(db=self._db, scorefile_paths=scorefile_paths)

    def import_match_tables(self, paths: Sequence[Pathish]) -> None:
        """Load the gated match table written by MATCH_COMBINE.

        `target_alleles=False` declares only the columns SCORE reads, which is
        what makes MATCH's ungated `matches.parquet` fail to load rather than
        score accessions that failed `min_overlap`.
        """
        stage.import_match_tables(db=self._db, paths=paths, target_alleles=False)

    def calculate_scores(self) -> None:
        with self._db.connect() as conn:
            create_score_table(conn=conn)

        sampleset = self._cache.sampleset
        # one shard per SCORE task, so ask for it before doing any work
        shard = self._cache.shard
        # the dosage array is (this shard's variants, this shard's samples), and the
        # same IDs label the scores calculated from it
        sample_ids = self._cache.sample_ids

        # the dosage array is the largest thing this stage creates and nothing
        # downstream reads it, so it goes in the scratch
        pgs_store = zarr.storage.LocalStore(self._scratch)
        pgs_group = zarr.open_group(store=pgs_store, mode="w", path=f"pgs/{sampleset}")

        # the match table was gated by MATCH_COMBINE, so every accession left
        # in it passed min_overlap and none of the rejected ones cost compute.
        # Scoped to this shard: an accession with no variants here contributes
        # nothing, and the samples it covers are scored by another shard
        # (every chromosome covers the same samples - check_sample_union)
        with self._db.read_only() as conn:
            accessions_to_calculate = [
                row[0]
                for row in conn.execute(
                    """
                    SELECT DISTINCT accession FROM allele_match_table
                    WHERE sampleset = ? AND zarr_group = ?
                      AND match_result.match_summary = 'matched'
                    ORDER BY accession
                    """,
                    [sampleset, shard.path],
                ).fetchall()
            ]

        if not accessions_to_calculate:
            logger.info(
                f"No matched variants in {shard.path!r}; exporting no contributions"
            )
            return

        # get wide effect weights (a matrix) for this shard's genotypes
        variant_metadata: pl.DataFrame = store_group_weight_arrays(
            db=self._db,
            sampleset=sampleset,
            pgs_group=pgs_group,
            accessions=accessions_to_calculate,
            zarr_group=shard.path,
        )

        # scoring only reads the database from here on, so one read-only
        # connection covers the rest of the shard
        with self._db.read_only() as conn:
            # calculate and adjust dosage with dask
            is_missing_zarr, dosage_zarr = make_dosage_arrays(
                genotypes_group=self._cache.root,
                df=variant_metadata,
                zarr_group=shard.path,
                pgs_store=pgs_store,
                sampleset=sampleset,
                n_samples=len(sample_ids),
                conn=conn,
                n_workers=self._threads,
                n_minimum_impute=self._n_minimum_samples_for_impute,
            )

            with dask.config.set(scheduler="threads", num_workers=self._threads):
                # calculate weighted sum of scores
                effect_weights: zarr.Array = cast(
                    "zarr.Array", pgs_group["weight_matrix"]
                )
                scores: npt.NDArray[np.float64] = calculate_scores(
                    dosage_array=dosage_zarr, effect_weights=effect_weights
                ).compute()

                _score_df = calculate_score_statistics(
                    conn=conn,
                    accessions=accessions_to_calculate,
                    dosage_array=dosage_zarr,
                    score=scores,
                    sampleset=sampleset,
                    sample_ids=sample_ids,
                    is_missing_array=is_missing_zarr,
                )

        # insert scores into duckdb for final cleanup and export
        with self._db.connect() as conn:
            insert_scores(conn=conn, _score_df=_score_df)
        logger.info(f"Score calculation {sampleset=} finished")

    def export_partial_scores(self, out_path: Pathish) -> None:
        """Write this shard's contribution for SCORE_COMBINE to sum.

        Columns are listed rather than starred to pin the artefact schema, and
        `score_avg` is left out because it is GENERATED and regenerated on the
        combined table. `zarr_group` is what lets SCORE_COMBINE re-assert the
        uniqueness constraint that parquet cannot carry.
        """
        logger.info(f"Exporting partial scores to {str(out_path)}")
        with self._db.connect() as conn:
            conn.execute(
                """
                COPY (
                    SELECT sampleset, accession, n_matched, sample_id,
                           allele_count, score, dosage_sum,
                           $zarr_group AS zarr_group
                    FROM score_table
                ) TO $path (FORMAT parquet, COMPRESSION zstd)
                """,
                {"zarr_group": self._cache.shard.path, "path": str(out_path)},
            )
