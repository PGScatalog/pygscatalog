"""SCORE_COMBINE: reduce every shard's partial scores into the final output.

Once SCORE is fanned out no task produces final output, so this is a single duckdb
reduction. It needs no awareness of how the cache was split: the reduction is
valid because every (sample, variant) pair appears in exactly one shard, which
MATCH's shape invariants enforce.
"""

from __future__ import annotations

import logging
import pathlib
from typing import TYPE_CHECKING

from rich import print

from pgscatalog.calc.cli._stage import (
    add_stage_args,
    check_paths_exist,
    setup_stage,
)
from pgscatalog.calc.lib.score.score_combine import ScoreCombineStage

if TYPE_CHECKING:
    import argparse

logger = logging.getLogger(__name__)


def add_combine_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--partial_scores",
        type=pathlib.Path,
        nargs="+",
        required=True,
        help="partial_scores.parquet from every shard's pgsc_calc score task",
    )
    add_stage_args(parser)
    parser.set_defaults(func=combine_cli)


def combine_cli(args: argparse.Namespace) -> None:
    out_dir = setup_stage(args)

    stage = ScoreCombineStage(
        max_memory_gb=args.max_memory_gb,
        threads=args.threads,
        out_dir=out_dir,
    )

    stage.combine_scores(check_paths_exist(args.partial_scores))
    stage.export_scores(out_path=out_dir)

    print("Finished calculating :tada: Goodbye!")
