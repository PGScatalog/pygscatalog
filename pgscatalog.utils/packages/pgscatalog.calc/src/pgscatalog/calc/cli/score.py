from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from rich.logging import RichHandler
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

from pgscatalog.calc import ScorePipeline

from .load import unzip_zarr

if TYPE_CHECKING:
    import argparse

logger = logging.getLogger(__name__)

logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)],
)


def score_cli(args: argparse.Namespace) -> None:
    # check the output directory is empty
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if any(files := [str(x) for x in args.out_dir.iterdir()]):
        logger.critical(f"{str(args.out_dir)} directory must be empty")
        raise FileExistsError(files)

    # check all zip files exist and unzip them
    for zip_file in args.zarr_zip_file:
        if not zip_file.exists():
            raise FileNotFoundError(f"{zip_file} does not exist")

        unzip_zarr(zip_file, args.out_dir)

    # zarr directory called genotypes.zarr now contains all variant and genotypes data

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        MofNCompleteColumn(),
        BarColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        refresh_per_second=2,
    ) as progress:
        stages = [
            "Load scoring files",
            "Match variants",
            "Export logs",
            "Calculate scores",
            "Export scores",
        ]
        tasks = {stage: progress.add_task(stage, total=1) for stage in stages}

        progress.print(":test_tube:", "pgsc_calc score", ":dna:", "starting...")

        pipeline = ScorePipeline(
            max_memory_gb=args.max_memory_gb,
            threads=args.threads,
            out_dir=args.out_dir,
        )

        progress.start_task(tasks["Load scoring files"])
        pipeline.load_scores(scorefile_paths=args.score_paths)
        progress.update(tasks["Load scoring files"], advance=1)

        progress.start_task(tasks["Match variants"])
        pipeline.match_variants(min_overlap=args.min_overlap)
        progress.update(tasks["Match variants"], advance=1)

        progress.start_task(tasks["Export logs"])
        pipeline.export_full_match_log(out_directory=args.out_dir / "logs")
        pipeline.export_summary_match_log(
            out_path=args.out_dir / "logs" / "summary.csv"
        )
        progress.update(tasks["Export logs"], advance=1)

        progress.start_task(tasks["Calculate scores"])
        pipeline.calculate_scores()
        progress.update(tasks["Calculate scores"], advance=1)

        progress.start_task(tasks["Export scores"])
        pipeline.export_scores(out_path=args.out_dir)
        progress.update(tasks["Export scores"], advance=1)

        progress.print("Finished calculating :tada: Goodbye!")
