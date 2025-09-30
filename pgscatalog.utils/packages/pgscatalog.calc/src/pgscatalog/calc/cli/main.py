from __future__ import annotations

import argparse
import pathlib

from pgscatalog.calc import GenomeFileType
from pgscatalog.calc.cli._utils import check_positive, zero_one_float
from pgscatalog.calc.cli.load import load_cli
from pgscatalog.calc.cli.score import score_cli


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="pgsc_calc",
        description="""
    Polygenic Score (PGS) Catalog Calculator. 
    
    A set of programs to apply PGS scoring files to new target genomes.
    """,
    )
    subparsers = parser.add_subparsers()

    parser_load = subparsers.add_parser(
        "load", help="Load genomes into a cache directory"
    )
    parser_load_args(parser_load)

    parser_score = subparsers.add_parser("score", help="Calculate polygenic scores")
    parser_score_args(parser_score)

    args = parser.parse_args()

    # dispatch args to correct function
    args.func(args)


def parser_score_args(parser: argparse.ArgumentParser) -> None:
    """Arguments for the score program"""
    parser.add_argument(
        "--cache_dir",
        type=pathlib.Path,
        dest="cache_dir",
        help="Genome cache directory.",
        required=True,
    )
    parser.add_argument(
        "--score_paths",
        type=pathlib.Path,
        nargs="+",
        help="PGS Catalog scoring file paths (processed with pgscatalog-format)",
        required=True,
    )
    parser.add_argument(
        "--out_dir",
        type=pathlib.Path,
        help="Directory path for output",
        required=True,
    )
    parser.add_argument(
        "--min_overlap",
        type=zero_one_float,
        dest="min_overlap",
        help="Minimum variant overlap",
        required=False,
        default=0.75,
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
    parser.set_defaults(func=score_cli)


def parser_load_args(parser: argparse.ArgumentParser) -> None:
    """Arguments for the load program"""
    parser.add_argument(
        "--cache_dir",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "pgscatalog-genome-cache",
        dest="cache_dir",
        help="A directory to store cached parquet files.",
    )
    parser.add_argument(
        "--target_genomes",
        type=str,
        nargs="+",
        help="Target genomes to load in the format PATH:CHROM",
        required=True,
    )
    parser.add_argument(
        "--sampleset",
        type=str,
        help="Human label for the target genomes (e.g. UKBiobank)",
        required=True,
    )
    parser.add_argument(
        "--format",
        type=GenomeFileType,
        choices=list(x.value for x in GenomeFileType),
        help="Target genome format",
        required=True,
    )
    parser.add_argument(
        "--score_paths",
        type=pathlib.Path,
        nargs="+",
        help="PGS Catalog scoring file paths",
        required=True,
    )
    parser.add_argument(
        "--update_cache", action="store_true", help="Update an existing cache"
    )
    parser.add_argument("--verbose", action="store_true", help="More logging!")
    parser.add_argument(
        "--workers",
        type=check_positive,
        default=1,
        help="Number of Python worker processes to use (adjust to be ~number of "
        "cores).",
    )
    parser.add_argument(
        "--batch_size",
        type=check_positive,
        default=10_000,
        help="Number of variants to buffer before writing to database",
    )
    parser.add_argument(
        "--timeout_seconds",
        type=check_positive,
        default=3600,
        help="Timeout for worker jobs in seconds. Helpful to kill workers that get "
        "stuck after crashing and not raising a valid Python exception",
    )
    parser.add_argument(
        "--bgen_sample_file",
        type=pathlib.Path,
        help="Path to a BGEN sample file",
        required=False,
    )
    parser.set_defaults(func=load_cli)


if __name__ == "__main__":
    main()
