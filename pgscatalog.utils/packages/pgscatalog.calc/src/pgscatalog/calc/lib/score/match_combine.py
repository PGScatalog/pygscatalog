"""MATCH_COMBINE: pool every chromosome's matches and apply `min_overlap`.

The match rate is genome-wide, so no single chromosome's MATCH task can decide it.
This stage has no cache to read - only the parquet artefacts MATCH published - so
the sampleset it is pooling has to be passed in.
"""

from __future__ import annotations

import logging
import pathlib
from typing import TYPE_CHECKING

from . import stage
from ._matchlog import add_match_log, get_ok_accessions
from .scoredb import ScoreDatabase

if TYPE_CHECKING:
    from collections.abc import Sequence

    import duckdb

    from pgscatalog.calc.lib.types import Pathish, PathishList


logger = logging.getLogger(__name__)


class MatchCombineStage:
    """The one MATCH_COMBINE task: every chromosome's match table, pooled.

    Parameters
    ----------
    out_dir : Pathish
        The task's output directory. `match_table.parquet` and `logs/` are published
        at its root; the database lives in its scratch.
    max_memory_gb : float
        Maximum RAM duckdb may use before intermediate results spill to disk.
    threads : int
        Maximum number of threads duckdb may use.
    sampleset : str
        The sampleset being pooled. There is no cache here to read it from.
    """

    def __init__(
        self,
        out_dir: Pathish,
        max_memory_gb: float,
        threads: int,
        sampleset: str,
    ):
        scratch = stage.scratch_dir(out_dir)
        self._sampleset = sampleset

        # every query this stage runs goes through the database, so the resource
        # limits the CLI was given reach all of them
        self._db = ScoreDatabase(
            path=scratch / "pgs.db",
            max_memory_gb=max_memory_gb,
            threads=threads,
        )

    def load_scores(self, scorefile_paths: Pathish | PathishList) -> None:
        """Copy scores processed with pgscatalog-format into a duckDB database table
        'scorevariants'
        """
        stage.load_scores(db=self._db, scorefile_paths=scorefile_paths)

    def import_match_tables(self, paths: Sequence[Pathish]) -> None:
        """Load the partial match tables written by every chromosome's MATCH task"""
        # the match log reports the target alleles, so this stage consumes them
        stage.import_match_tables(db=self._db, paths=paths, target_alleles=True)

    def combine_matches(self, min_overlap: float, out_path: Pathish) -> None:
        """Pool every chromosome's matches, apply min_overlap, and write the
        genome-wide match table restricted to accessions that passed"""
        if not (0 < min_overlap <= 1):
            raise ValueError(
                f"min_overlap must be positive and less than 1, got {min_overlap}"
            )

        with self._db.connect() as conn:
            add_match_log(
                conn=conn,
                min_overlap=min_overlap,
                sampleset=self._sampleset,
            )

        with self._db.read_only() as conn:
            accessions = get_ok_accessions(conn=conn, sampleset=self._sampleset)

        logger.info(f"Exporting gated match table to {str(out_path)}")
        with self._db.connect() as conn:
            conn.execute(
                """
                COPY (
                    SELECT * EXCLUDE (target_ref, target_alts)
                    FROM allele_match_table
                    WHERE accession IN $accessions
                ) TO $path (FORMAT parquet, COMPRESSION zstd)
                """,
                {"accessions": accessions, "path": str(out_path)},
            )

    def export_full_match_log(self, out_directory: Pathish) -> None:
        logger.info(f"Exporting full match log to {str(out_directory)}")
        pathlib.Path(out_directory).mkdir(parents=True, exist_ok=True)

        with self._db.connect() as conn:
            conn.table("score_log_table").order("accession, row_nr").write_csv(
                file_name=str(out_directory),
                partition_by=["sampleset", "accession", "chr_name"],
                write_partition_columns=True,
                compression="gzip",
                overwrite=True,
            )
        logger.info("Full match log exported")

    def export_summary_match_log(self, out_path: Pathish) -> None:
        logger.info(f"Exporting summary match log to {str(out_path)}")
        pathlib.Path(out_path).parent.mkdir(parents=True, exist_ok=True)

        with self._db.connect() as conn:
            (
                conn.table("summary_log_table")
                .order("sampleset,accession,match_summary")
                .write_csv(file_name=str(out_path))
            )

        logger.info("Summary match log exported")


def check_sample_union(
    conn: duckdb.DuckDBPyConnection, paths: Sequence[Pathish]
) -> None:
    """Every chromosome in the cache must cover the same samples, except Y.

    Y is the only chromosome a person can lack entirely - everyone has at least one
    X and everyone has MT - so a male-only chrY leaf covers a subset of the sampleset
    and is legitimate rather than truncated (ADR 0007). The two are indistinguishable,
    which is the known cost of admitting the first.

    Scoped to chromosomes present in the cache, which is automatic: a chromosome
    whose index query returned nothing is deleted by `cleanup_empty_zarr_group`,
    runs no MATCH task and writes no samples file.
    """
    mismatched = conn.execute(
        """
        WITH samples AS (SELECT * FROM read_parquet($paths)),
             union_of_all AS (SELECT DISTINCT sampleset, sample_id FROM samples),
             chromosomes AS (SELECT DISTINCT sampleset, chromosome FROM samples
                             WHERE chromosome != 'Y'),
             expected AS (SELECT * FROM chromosomes JOIN union_of_all
                          USING (sampleset))
        SELECT DISTINCT chromosome
        FROM (SELECT * FROM expected EXCEPT SELECT * FROM samples)
        """,
        {"paths": [str(p) for p in paths]},
    ).fetchall()

    if mismatched:
        raise ValueError(
            f"Chromosomes {sorted(x[0] for x in mismatched)} are missing samples "
            "that other chromosomes have. Every chromosome of a sampleset must "
            "cover the same samples, or scores would be summed over a different "
            "cohort per chromosome."
        )
