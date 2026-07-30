"""Shared plumbing for the four stages of split-apply-combine scoring.

MATCH -> MATCH_COMBINE -> SCORE -> SCORE_COMBINE each run as their own Nextflow
task (ADR 0001), so each is its own subcommand with the same output directory
and resource options.
"""

from __future__ import annotations

import logging
import pathlib
from typing import TYPE_CHECKING

from rich.logging import RichHandler

from pgscatalog.calc.cli._utils import check_positive

if TYPE_CHECKING:
    import argparse

logger = logging.getLogger(__name__)


def add_stage_args(parser: argparse.ArgumentParser) -> None:
    """Options every stage takes"""
    parser.add_argument(
        "--out_dir",
        type=pathlib.Path,
        help="Directory path for output",
        required=True,
    )
    parser.add_argument(
        "--threads",
        type=check_positive,
        help="Number of threads to use",
        required=False,
        default=1,
    )
    parser.add_argument(
        "--max_memory_gb",
        # have to double %: argparse uses % to parse strings
        help="Maximum RAM to use in GB (default: ~90%% of system RAM if not set)",
        type=check_positive,
        required=False,
        default=16,
    )
    parser.add_argument("--verbose", action="store_true", help="More logging!")


def add_score_paths_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--score_paths",
        type=pathlib.Path,
        nargs="+",
        help="PGS Catalog scoring file paths (processed with pgscatalog-format)",
        required=True,
    )


def setup_stage(args: argparse.Namespace) -> pathlib.Path:
    """Configure logging and prepare an empty output directory"""
    logging.basicConfig(
        level="INFO" if args.verbose else "WARNING",
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True)],
    )

    out_dir: pathlib.Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # the output directory contains results, so don't overwrite stuff
    if files := [str(x) for x in out_dir.iterdir()]:
        raise FileExistsError(f"{str(out_dir)} directory must be empty: {files}")

    return out_dir


def check_paths_exist(paths: list[pathlib.Path]) -> list[pathlib.Path]:
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(f"{path} does not exist")
    return paths
