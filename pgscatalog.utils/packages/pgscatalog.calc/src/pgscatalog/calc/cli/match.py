"""MATCH: match scoring file variants against one chromosome's shards.

Fanning out per chromosome rather than per shard is the load-bearing choice: a
chromosome-level task holds every shard's variant metadata for that chromosome,
so it matches against a complete variant list and no shard has a partial view. It
reads 1D metadata only, never the 3D genotype arrays, so peak memory doesn't
scale with sample count.
"""

from __future__ import annotations

import logging
import pathlib
from typing import TYPE_CHECKING

from rich import print

from pgscatalog.calc.cli._stage import (
    add_score_paths_arg,
    add_stage_args,
    check_paths_exist,
    setup_stage,
)
from pgscatalog.calc.cli._utils import check_positive
from pgscatalog.calc.lib.config import config
from pgscatalog.calc.lib.score.match import MatchStage
from pgscatalog.calc.lib.score.stage import scratch_dir

from .load import unzip_zarr

if TYPE_CHECKING:
    import argparse

logger = logging.getLogger(__name__)


def add_match_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--zarr_zip_file",
        type=pathlib.Path,
        nargs="+",
        dest="zarr_zip_file",
        required=True,
        help="Zarr zip files created by pgsc_calc load: every shard of one "
        "chromosome, and only that chromosome.",
    )
    add_score_paths_arg(parser)
    parser.add_argument(
        "--keep_multiallelic",
        dest="keep_multiallelic",
        action="store_true",
        help="Flag to allow matching to multiallelic variants (default: false).",
        required=False,
    )
    parser.add_argument(
        "--keep_ambiguous",
        dest="keep_ambiguous",
        action="store_true",
        help="""
            Flag to force the program to keep variants with
            ambiguous alleles, (e.g. A/T and G/C SNPs), which are normally
            excluded (default: false). In this case the program proceeds
            assuming that the genotype data is on the same strand as the
            GWAS whose summary statistics were used to construct the score.
        """,
    )
    parser.add_argument(
        "--ignore_strand_flips",
        dest="match_flipped",
        action="store_false",
        help="""
            Flag to not consider matched variants that may be reported on the
            opposite strand. Default behaviour is to check whether unmatched
            variants match a reverse complemented reading of themselves. Variants
            that only match that way are excluded, and the match log says so.
        """,
        required=False,
    )
    parser.add_argument(
        "--batch_size",
        type=check_positive,
        default=5_000,
        help="Number of variants to buffer before writing to zarr",
    )
    add_stage_args(parser)
    parser.set_defaults(func=match_cli)


def match_cli(args: argparse.Namespace) -> None:
    out_dir = setup_stage(args)
    config.ZARR_VARIANT_CHUNK_SIZE = args.batch_size

    # the cache is scratch: this task publishes two parquet files and nothing else
    for zip_file in check_paths_exist(args.zarr_zip_file):
        unzip_zarr(zip_file, scratch_dir(out_dir))

    stage = MatchStage(
        max_memory_gb=args.max_memory_gb,
        threads=args.threads,
        out_dir=out_dir,
    )

    stage.load_scores(scorefile_paths=args.score_paths)

    # no min_overlap here: the match rate is genome-wide, so only MATCH_COMBINE
    # can decide it
    chromosome = stage.match_variants(
        match_multiallelic=args.keep_multiallelic,
        match_ambiguous=args.keep_ambiguous,
        match_flipped=args.match_flipped,
    )

    stage.export_match_table(out_dir / "matches.parquet")
    stage.export_samples(out_dir / "samples.parquet", chromosome=chromosome)

    print(f"Finished matching chromosome {chromosome} :tada:")
