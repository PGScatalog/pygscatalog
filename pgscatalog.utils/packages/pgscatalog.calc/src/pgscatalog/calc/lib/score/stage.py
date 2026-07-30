"""What the four stages of split-apply-combine share.

MATCH, MATCH_COMBINE, SCORE and SCORE_COMBINE are four modules, one per stage, each
taking only what it needs over the cache seam (docs/tickets/0012) and the database
seam (docs/tickets/0016). They are not variations of one thing, so there is no base
class here: what genuinely repeats is the directory layout, and loading the two
tables more than one stage loads.

A stage's output directory holds two different kinds of thing::

    out_dir/
      matches.parquet     <- published: a Nextflow process declares it by name
      samples.parquet
      scratch/            <- everything the stage only needs while it runs
        pgs.db
        genotypes.zarr/
        pgs/{sampleset}/dosage

Published artefacts stay at the root because the process declarations and
`publishDir` globs in `pgsc_calc` name them there (ADR 0001). Everything else -
including the dosage array, the largest thing SCORE creates - goes under `scratch`,
so a task's outputs are enumerable without listing its scratch (docs/tickets/0017).
"""

from __future__ import annotations

import logging
import pathlib
from typing import TYPE_CHECKING

from pgscatalog.calc.lib.scorefile import load_scoring_files

from ._matchvariants import _setup_enums, create_allele_match_table

if TYPE_CHECKING:
    from collections.abc import Sequence

    from pgscatalog.calc.lib.types import Pathish, PathishList

    from .scoredb import ScoreDatabase

logger = logging.getLogger(__name__)


def scratch_dir(out_dir: Pathish) -> pathlib.Path:
    """Where a stage keeps what it isn't publishing, created if it doesn't exist.

    Parameters
    ----------
    out_dir : Pathish
        The stage's output directory, which must already exist: the published
        artefacts land at its root.
    """
    directory = pathlib.Path(out_dir).resolve()

    # the output directory contains results, so don't overwrite stuff
    if not directory.is_dir():
        raise NotADirectoryError(f"{str(directory)} must be a directory")

    scratch = directory / "scratch"
    scratch.mkdir(exist_ok=True)
    return scratch


def load_scores(db: ScoreDatabase, scorefile_paths: Pathish | PathishList) -> None:
    """Copy scores processed with pgscatalog-format into a duckDB database table
    'scorevariants'

    Three of the four stages do this: MATCH and SCORE to match and weight variants,
    MATCH_COMBINE because the `min_overlap` denominator is every variant in the
    scoring file, including variants on chromosomes that aren't in the cache at all.

    Parameters
    ----------
    db : ScoreDatabase
        The stage's score database.
    scorefile_paths : Pathish | PathishList
        One scoring file path, or a list of them.
    """
    if not isinstance(scorefile_paths, list):
        scorefile_paths = [scorefile_paths]

    with db.connect() as conn:
        load_scoring_files(conn=conn, scorefile_paths=scorefile_paths)


def import_match_tables(db: ScoreDatabase, paths: Sequence[Pathish]) -> None:
    """Load partial match tables written by MATCH (or the gated table written
    by MATCH_COMBINE, which drops the target allele columns)

    MATCH_COMBINE pools every chromosome's partial table; SCORE reads back the one
    gated table. Neither of them runs the matcher, so the enums the match table
    declares have to be set up here.

    Parameters
    ----------
    db : ScoreDatabase
        The stage's score database.
    paths : Sequence[Pathish]
        Parquet match tables to load.
    """
    logger.info(f"Importing {len(paths)} match tables")
    with db.connect() as conn:
        # score_log_table declares enum-typed columns, and update_match_table
        # (the only other caller) never runs in a combine task
        _setup_enums(conn)
        create_allele_match_table(conn)
        conn.execute(
            """
            INSERT INTO allele_match_table BY NAME
            SELECT * FROM read_parquet($paths)
            """,
            {"paths": [str(p) for p in paths]},
        )
