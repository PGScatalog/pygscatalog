"""Every connection to the score database carries the stage's resource limits
(docs/tickets/0016).

duckdb's default is a share of *visible* RAM, which inside a scheduled task is
usually the host's rather than the task's limit, so a query that ignores
`--max_memory_gb` gets the task killed by the scheduler instead of spilling to disk.
Each stage is a separate Nextflow task with its own limits (ADR 0001), so this is
asserted by running the four stages and reading the settings back from inside them,
rather than by grepping for connection calls.
"""

import duckdb

from pgscatalog.calc.lib.score.scoredb import ScoreDatabase

from .conftest import run_pipeline

MAX_MEMORY_GB = "3"
THREADS = "2"


def _settings(conn: duckdb.DuckDBPyConnection) -> tuple:
    return conn.execute(
        "SELECT current_setting('max_memory'), current_setting('threads')"
    ).fetchone()


def test_every_connection_carries_the_resource_policy(
    grid_targets, scorefile, tmp_path, monkeypatch
):
    """MATCH, MATCH_COMBINE, SCORE and SCORE_COMBINE all hand out connections from
    ScoreDatabase, so every query they run is bounded by what the CLI was given -
    including the PIVOT in create_wide_weights_table, the memory-hungriest query in
    SCORE and the site of the alpha1 out-of-memory failure (pgsc_calc#459).
    """
    observed: list[tuple] = []
    open_connection = ScoreDatabase._connect

    def spy(self, *, read_only):
        conn = open_connection(self, read_only=read_only)
        observed.append(_settings(conn))
        return conn

    monkeypatch.setattr(ScoreDatabase, "_connect", spy)

    run_pipeline(
        grid_targets,
        tmp_path,
        scorefile,
        monkeypatch,
        max_memory_gb=MAX_MEMORY_GB,
        threads=THREADS,
    )

    with duckdb.connect(
        config={"max_memory": f"{MAX_MEMORY_GB}GB", "threads": THREADS}
    ) as reference:
        expected = _settings(reference)

    # a four-shard run opens a connection per phase of every stage, so anything much
    # smaller than this means a stage stopped going through ScoreDatabase
    assert len(observed) > 20
    assert set(observed) == {expected}
