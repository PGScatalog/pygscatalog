"""SCORE: calculate partial polygenic scores for one shard.

Scoring used to merge every shard's genotypes into a single zarr store and compute
over it as one dask array in one task, which is the OOM in pgsc_calc#459. Fanning
out per shard is the fix and the point of alpha2: peak memory scales with one
shard rather than the sampleset.

Accessions rejected by `min_overlap` were dropped by MATCH_COMBINE, so they cost
no scoring compute here.
"""

from __future__ import annotations

import logging
import pathlib
from typing import TYPE_CHECKING

from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

from pgscatalog.calc.cli._stage import (
    add_score_paths_arg,
    add_stage_args,
    check_paths_exist,
    setup_stage,
)
from pgscatalog.calc.cli._utils import check_positive
from pgscatalog.calc.lib.config import config
from pgscatalog.calc.lib.score.score import ScoreStage
from pgscatalog.calc.lib.score.stage import scratch_dir

from .load import unzip_zarr

if TYPE_CHECKING:
    import argparse

logger = logging.getLogger(__name__)


def add_score_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--zarr_zip_file",
        type=pathlib.Path,
        required=True,
        help="A zarr zip file created by pgsc_calc load, holding one shard.",
    )
    add_score_paths_arg(parser)
    parser.add_argument(
        "--match_table",
        type=pathlib.Path,
        required=True,
        help="match_table.parquet written by pgsc_calc match-combine",
    )
    parser.add_argument(
        "--batch_size",
        type=check_positive,
        default=5_000,
        help="Number of variants to buffer before writing to zarr",
    )
    add_stage_args(parser)
    parser.set_defaults(func=score_cli)


def score_cli(args: argparse.Namespace) -> None:
    out_dir: pathlib.Path = setup_stage(args)

    logger.info(f"Default value {config.ZARR_VARIANT_CHUNK_SIZE=}")
    config.ZARR_VARIANT_CHUNK_SIZE = args.batch_size
    logger.info(f"After update {config.ZARR_VARIANT_CHUNK_SIZE=}")

    # the cache is scratch, and so are the dosage arrays derived from it: this task
    # publishes partial_scores.parquet and nothing else
    (zip_file,) = check_paths_exist([args.zarr_zip_file])
    unzip_zarr(zip_file, scratch_dir(out_dir))

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        MofNCompleteColumn(),
        BarColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        refresh_per_second=2,
    ) as progress:
        phases = [
            "Load scoring files",
            "Load match table",
            "Calculate scores",
            "Export partial scores",
        ]
        tasks = {phase: progress.add_task(phase, total=1) for phase in phases}

        progress.print(":test_tube:", "pgsc_calc score", ":dna:", "starting...")

        stage = ScoreStage(
            max_memory_gb=args.max_memory_gb,
            threads=args.threads,
            out_dir=out_dir,
        )

        progress.start_task(tasks["Load scoring files"])
        stage.load_scores(scorefile_paths=args.score_paths)
        progress.update(tasks["Load scoring files"], advance=1)

        progress.start_task(tasks["Load match table"])
        stage.import_match_tables(check_paths_exist([args.match_table]))
        progress.update(tasks["Load match table"], advance=1)

        progress.start_task(tasks["Calculate scores"])
        stage.calculate_scores()
        progress.update(tasks["Calculate scores"], advance=1)

        progress.start_task(tasks["Export partial scores"])
        stage.export_partial_scores(out_dir / "partial_scores.parquet")
        progress.update(tasks["Export partial scores"], advance=1)

        progress.print("Finished scoring this shard :tada:")
