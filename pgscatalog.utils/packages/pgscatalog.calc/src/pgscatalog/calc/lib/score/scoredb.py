"""The score database, and the resources its queries are allowed to spend.

Every stage of split-apply-combine works over one duckdb database (`pgs.db`) in the
stage's output directory. Each stage is a separate Nextflow task with its own
resource limits (ADR 0001), so the limits the CLI was given have to reach every
query that stage makes: duckdb's default is a share of *visible* RAM, which inside
a scheduled task is usually the host's rather than the task's, so exceeding it gets
the task killed by the scheduler instead of spilled to disk.

This module is the only place that opens a connection to it. Functions that query
the database take a connection, never a path (docs/tickets/0016).
"""

from __future__ import annotations

import logging
import pathlib
from typing import TYPE_CHECKING

import duckdb

if TYPE_CHECKING:
    from pgscatalog.calc.lib.types import Pathish

logger = logging.getLogger(__name__)


class ScoreDatabase:
    """One stage's score database: where it lives and what it may spend.

    Parameters
    ----------
    path : Pathish
        Path to the duckdb database file (`pgs.db` in the stage's output directory).
    max_memory_gb : float
        Maximum RAM duckdb may use before intermediate results spill to disk.
    threads : int
        Maximum number of threads duckdb may use.
    """

    def __init__(self, path: Pathish, max_memory_gb: float, threads: int):
        self._path = pathlib.Path(path)
        # note: duckdb will use 90% of available memory and all threads by default
        self._config: dict[str, str] = {
            "max_memory": f"{max_memory_gb}GB",
            "threads": str(threads),
        }

    @property
    def path(self) -> pathlib.Path:
        return self._path

    def connect(self) -> duckdb.DuckDBPyConnection:
        """Open a writable connection carrying this stage's resource limits"""
        return self._connect(read_only=False)

    def read_only(self) -> duckdb.DuckDBPyConnection:
        """Open a read-only connection carrying this stage's resource limits.

        Queries that only read take one of these, so a query that grows a write by
        accident fails instead of quietly editing a database that holds results.
        """
        return self._connect(read_only=True)

    def _connect(self, *, read_only: bool) -> duckdb.DuckDBPyConnection:
        # ponytail: connections are opened and closed around each phase rather than
        # pooled. duckdb keeps one instance per database file per process and refuses
        # a second connection asking for a different configuration, so a read-only
        # connection can only be opened while no writable one is held open.
        logger.info(f"Connecting to {self._path} {read_only=} {self._config}")
        return duckdb.connect(str(self._path), read_only=read_only, config=self._config)
