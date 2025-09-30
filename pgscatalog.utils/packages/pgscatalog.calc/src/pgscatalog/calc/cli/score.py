from __future__ import annotations

import argparse
import logging
import pathlib

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

logger = logging.getLogger("pgscatalog.wgs")

logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)],
)


def score_cli(args: argparse.Namespace) -> None:
    out_dir = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

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

        progress.print(
            ":test_tube:", "pgscatalog-calculate-scores", ":dna:", "starting..."
        )

        pipeline = ScorePipeline(
            cache_dir=args.cache_dir,
            max_memory_gb=args.max_memory_gb,
            threads=args.threads,
        )

        progress.start_task(tasks["Load scoring files"])
        pipeline.load_scores(scorefile_paths=args.score_paths)
        progress.update(tasks["Load scoring files"], advance=1)

        progress.start_task(tasks["Match variants"])
        pipeline.match_variants(min_overlap=args.min_overlap)
        progress.update(tasks["Match variants"], advance=1)

        progress.start_task(tasks["Export logs"])
        pipeline.export_full_match_log(out_directory=out_dir / "logs")
        pipeline.export_summary_match_log(out_path=out_dir / "logs" / "summary.csv")
        progress.update(tasks["Export logs"], advance=1)

        progress.start_task(tasks["Calculate scores"])
        pipeline.calculate_scores()
        progress.update(tasks["Calculate scores"], advance=1)

        progress.start_task(tasks["Export scores"])
        pipeline.export_scores(out_path=out_dir / "scores")
        progress.update(tasks["Export scores"], advance=1)

        progress.print("Finished calculating :tada: Goodbye!")
