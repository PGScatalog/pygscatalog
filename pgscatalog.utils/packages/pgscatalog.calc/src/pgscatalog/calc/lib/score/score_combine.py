"""SCORE_COMBINE: sum every shard's partial scores into the final output.

The reduction is valid without any awareness of how the cache was split, because
every (sample, variant) pair appears in exactly one shard (docs/tickets/0008). So
this stage has no cache, no shard, no sampleset and no scoring files: parquet in,
`scores.txt.gz` out.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from . import stage
from ._pgs import create_score_table, write_scores
from .scoredb import ScoreDatabase

if TYPE_CHECKING:
    from collections.abc import Sequence

    import duckdb

    from pgscatalog.calc.lib.types import Pathish


logger = logging.getLogger(__name__)


class ScoreCombineStage:
    """The one SCORE_COMBINE task: every shard's partial scores, summed.

    Parameters
    ----------
    out_dir : Pathish
        The task's output directory. `scores.txt.gz` is published at its root; the
        database lives in its scratch.
    max_memory_gb : float
        Maximum RAM duckdb may use before intermediate results spill to disk.
    threads : int
        Maximum number of threads duckdb may use.
    """

    def __init__(self, out_dir: Pathish, max_memory_gb: float, threads: int):
        scratch = stage.scratch_dir(out_dir)

        # every query this stage runs goes through the database, so the resource
        # limits the CLI was given reach all of them (docs/tickets/0016)
        self._db = ScoreDatabase(
            path=scratch / "pgs.db",
            max_memory_gb=max_memory_gb,
            threads=threads,
        )

    def combine_scores(self, paths: Sequence[Pathish]) -> None:
        """Sum every shard's partial scores.

        Valid without any awareness of how the cache was split, because every
        (sample, variant) pair appears in exactly one shard (docs/tickets/0008).
        """
        logger.info(f"Combining {len(paths)} partial score tables")
        parquet_paths = [str(p) for p in paths]

        with self._db.connect() as conn:
            create_score_table(conn=conn)

            # n_rows / n_shards ride along on the aggregate that has to happen
            # anyway, so checking costs a pass over the combined table rather than
            # a second pass over the partials
            conn.execute(
                """
                CREATE OR REPLACE TEMP TABLE combined AS
                SELECT sampleset,
                       accession,
                       sample_id,
                       -- already genome-wide (COUNT(DISTINCT row_nr) over the
                       -- pooled match table), so carried through, never summed
                       any_value(n_matched)       AS n_matched,
                       SUM(allele_count)          AS allele_count,
                       SUM(score)                 AS score,
                       SUM(dosage_sum)            AS dosage_sum,
                       COUNT(*)                   AS n_rows,
                       COUNT(DISTINCT zarr_group) AS n_shards
                FROM read_parquet($paths)
                GROUP BY sampleset, accession, sample_id
                """,
                {"paths": parquet_paths},
            )

            check_shards_contribute_once(conn=conn, paths=parquet_paths)

            # column order follows the score_table DDL; score_avg is generated
            conn.execute("""
                INSERT INTO score_table
                SELECT sampleset, accession, n_matched, sample_id,
                       allele_count, score, dosage_sum
                FROM combined
                """)

    def export_scores(self, out_path: Pathish) -> None:
        with self._db.connect() as conn:
            write_scores(conn=conn, out_dir=out_path)


def check_shards_contribute_once(
    conn: duckdb.DuckDBPyConnection, paths: Sequence[Pathish]
) -> None:
    """Each shard must reach SCORE_COMBINE exactly once.

    Reads the `combined` temp table built by `combine_scores`, where a sample's
    partial scores have already been counted: they must come from as many distinct
    shards as there are rows to sum.

    This is the constraint the `score_table` primary key cannot express. Within a
    shard the key already guarantees one row per (accession, sample), but a
    duplicate arrives as a separate file - a glob collecting the same output twice,
    or a retried task leaving its earlier output alongside the new one - and parquet
    carries no keys across that boundary. Summing it twice doubles that shard's
    score, allele count and dosage sum with no error and no NaN.
    """
    duplicated = conn.execute("""
        SELECT accession, sample_id, n_rows, n_shards
        FROM combined
        WHERE n_rows <> n_shards
        ORDER BY accession, sample_id
        LIMIT 3
        """).fetchall()

    if not duplicated:
        return

    # only on the error path, so the happy path pays nothing to name the shard
    accession, sample_id, _, _ = duplicated[0]
    shards = conn.execute(
        """
        SELECT zarr_group, COUNT(*) AS n
        FROM read_parquet($paths)
        WHERE accession = $accession AND sample_id = $sample_id
        GROUP BY zarr_group
        HAVING COUNT(*) > 1
        ORDER BY zarr_group
        """,
        {
            "paths": [str(p) for p in paths],
            "accession": accession,
            "sample_id": sample_id,
        },
    ).fetchall()

    listed = ", ".join(f"{group} ({n} times)" for group, n in shards)
    affected = ", ".join(f"{acc}/{sample}" for acc, sample, _, _ in duplicated)
    raise ValueError(
        f"Shard {listed} reached SCORE_COMBINE more than once. Partial scores are "
        "summed per sample, so scoring a shard twice doubles its score, allele "
        f"count and dosage sum. Affected (accession/sample): {affected}. Check that "
        "each shard contributes exactly one partial_scores.parquet - a stale output "
        "from a retried task, or a glob matching the same file twice, both do this."
    )
