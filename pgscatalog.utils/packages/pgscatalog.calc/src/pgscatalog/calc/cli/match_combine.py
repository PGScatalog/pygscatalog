"""MATCH_COMBINE: pool every chromosome's matches and decide which accessions
pass `min_overlap`.

No single chromosome's MATCH task can make this decision, because the match rate
is genome-wide (docs/tickets/0004). Pooling happens once, centrally, and its
output is the gated match table that SCORE fans out over.
"""

from __future__ import annotations

import logging
import pathlib
from typing import TYPE_CHECKING

import duckdb
from rich import print

from pgscatalog.calc.cli._stage import (
    add_score_paths_arg,
    add_stage_args,
    check_paths_exist,
    setup_stage,
)
from pgscatalog.calc.cli._utils import zero_one_float
from pgscatalog.calc.lib.score.match_combine import (
    MatchCombineStage,
    check_sample_union,
)

if TYPE_CHECKING:
    import argparse

logger = logging.getLogger(__name__)


def add_match_combine_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--matches",
        type=pathlib.Path,
        nargs="+",
        required=True,
        help="matches.parquet from every chromosome's pgsc_calc match task",
    )
    parser.add_argument(
        "--samples",
        type=pathlib.Path,
        nargs="+",
        required=True,
        help="samples.parquet from every chromosome's pgsc_calc match task",
    )
    add_score_paths_arg(parser)
    parser.add_argument(
        "--sampleset",
        type=str,
        required=True,
        help="Human label for the target genomes (e.g. UKBiobank)",
    )
    parser.add_argument(
        "--min_overlap",
        type=zero_one_float,
        dest="min_overlap",
        help="Minimum variant overlap",
        required=False,
        default=0.75,
    )
    add_stage_args(parser)
    parser.set_defaults(func=match_combine_cli)


def match_combine_cli(args: argparse.Namespace) -> None:
    out_dir = setup_stage(args)

    # there is no cache here to read the sampleset from, only the parquet files
    # MATCH published, so it's passed in
    stage = MatchCombineStage(
        max_memory_gb=args.max_memory_gb,
        threads=args.threads,
        out_dir=out_dir,
        sampleset=args.sampleset,
    )

    with duckdb.connect() as conn:
        check_sample_union(conn=conn, paths=check_paths_exist(args.samples))

    # reloading the scoring files is required, not redundant: the min_overlap
    # denominator is every variant in the scoring file, including variants on
    # chromosomes that aren't in the cache at all
    stage.load_scores(scorefile_paths=args.score_paths)
    stage.import_match_tables(check_paths_exist(args.matches))

    stage.combine_matches(
        min_overlap=args.min_overlap, out_path=out_dir / "match_table.parquet"
    )
    stage.export_full_match_log(out_directory=out_dir / "logs")
    stage.export_summary_match_log(out_path=out_dir / "logs" / "summary.csv")

    print("Finished combining matches :tada:")
