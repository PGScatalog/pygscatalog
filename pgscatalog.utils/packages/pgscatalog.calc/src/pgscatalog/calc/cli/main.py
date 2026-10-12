from __future__ import annotations

import argparse
import pathlib

from pgscatalog.core.cli.download_cli import add_download_args, download_cli
from pgscatalog.core.cli.format_cli import add_format_args, format_cli

from pgscatalog.calc.cli._utils import check_positive
from pgscatalog.calc.cli.combine import add_combine_args
from pgscatalog.calc.cli.load import add_load_input_args, load_cli
from pgscatalog.calc.cli.match import add_match_args
from pgscatalog.calc.cli.match_combine import add_match_combine_args
from pgscatalog.calc.cli.run import add_run_args
from pgscatalog.calc.cli.score import add_score_args


def build_parser() -> argparse.ArgumentParser:
    """The parser for every subcommand

    A serial run composes argv and parses it here rather than calling stage
    functions directly, so argparse checks that `run` still speaks each stage's
    documented argument language (ADR 0013).
    """
    parser = argparse.ArgumentParser(
        prog="pgsc_calc",
        description="""
    Polygenic Score (PGS) Catalog Calculator.

    A set of programs to apply PGS scoring files to new target genomes.

    Scoring is split into four stages which fan out at different levels of the
    genotype cache: match (per chromosome) -> match-combine (once) ->
    score (per shard) -> combine (once).
    """,
    )
    subparsers = parser.add_subparsers(required=True)

    parser_load = subparsers.add_parser(
        "load", help="Load genomes into a cache directory"
    )
    parser_load_args(parser_load)

    add_match_args(
        subparsers.add_parser(
            "match", help="Match scoring file variants against one chromosome"
        )
    )
    add_match_combine_args(
        subparsers.add_parser(
            "match-combine",
            help="Pool match counts across chromosomes and apply min_overlap",
        )
    )
    add_score_args(
        subparsers.add_parser("score", help="Calculate polygenic scores for one shard")
    )
    add_combine_args(
        subparsers.add_parser("combine", help="Sum partial scores from every shard")
    )
    add_run_args(
        subparsers.add_parser(
            "run", help="Run every scoring stage in sequence, on one machine"
        )
    )

    parser_format = subparsers.add_parser(
        "format", help="Format scoring files to a consistent structure "
    )
    parser_format_args(parser_format)

    parser_download = subparsers.add_parser(
        "download", help="Download PGS Catalog scoring files "
    )
    parser_download_args(parser_download)

    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)

    # dispatch args to correct function
    args.func(args)


def parser_download_args(parser: argparse.ArgumentParser) -> None:
    add_download_args(parser)
    parser.set_defaults(func=download_cli)


def parser_format_args(parser: argparse.ArgumentParser) -> None:
    add_format_args(parser)
    parser.set_defaults(func=format_cli)


def parser_load_args(parser: argparse.ArgumentParser) -> None:
    """Arguments for the load program"""
    parser.add_argument(
        "--cache_dir",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "pgscatalog-genome-cache",
        dest="cache_dir",
        help="A directory to store the cache.",
    )
    parser.add_argument(
        "--zarr_zip_file",
        type=pathlib.Path,
        dest="zarr_zip_file",
        required=False,
        help="An existing zarr cache to update.",
    )
    add_load_input_args(parser)
    parser.add_argument("--verbose", action="store_true", help="More logging!")
    parser.add_argument(
        "--threads",
        type=check_positive,
        default=1,
        help="Number of dask threads to use (between 2 - 4 is good)",
    )
    parser.set_defaults(func=load_cli)


if __name__ == "__main__":
    main()
