import argparse
import concurrent.futures
import json
import logging
import pathlib
import sys
import textwrap

from pgscatalog.corelib import GenomeBuild, ScoringFile

from pgscatalog.combineapp._combine import normalise, get_variant_log, TextFileWriter

logger = logging.getLogger(__name__)


def run():
    logging.basicConfig(
        format="%(asctime)s %(name)s %(levelname)-8s %(message)s",
        level=logging.WARNING,
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    args = parse_args()

    if args.verbose:
        logger.setLevel(logging.INFO)

    out_path = pathlib.Path(args.outfile)

    if out_path.exists():
        raise FileExistsError(f"{args.outfile}")

    match x := out_path.name:
        case _ if x.endswith("gz"):
            compress_output = True
        case _:
            compress_output = False

    paths = list(set(args.scorefiles))  # unique paths only
    scoring_files = [ScoringFile(x) for x in paths]
    target_build = GenomeBuild.from_string(args.target_build)

    for x in scoring_files:
        if x.genome_build != target_build:
            raise NotImplementedError(
                f"{x.pgs_id} build {x.genome_build} doesn't match {target_build}. "
                f"Liftover not implemented yet"
            )

    bad_builds = [x.pgs_id for x in scoring_files if x.genome_build != target_build]

    for bad_file in bad_builds:
        raise ValueError(f"{bad_file} doesn't match {target_build}, can't combine")
    else:
        logger.info(f"All builds match target build {target_build}")

    variant_log = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []
        for scorefile in scoring_files:
            logger.info(f"Submitting {scorefile!r}")
            futures.append(
                executor.submit(
                    normalise,
                    scorefile,
                )
            )

        for future in concurrent.futures.as_completed(futures):
            writer = TextFileWriter(compress=compress_output, filename=out_path)
            normalised_score = future.result()
            writer.write(normalised_score)
            variant_log.append(get_variant_log(normalised_score))

    # TODO: fix drop_missing argument
    score_log = []
    for sf, log in zip(scoring_files, variant_log, strict=True):
        score_log.append(sf.get_log(variant_log=log))

    log_out_path = pathlib.Path(args.outfile).parent / args.logfile
    with open(log_out_path, "w") as f:
        logger.info(f"Writing log to {f.name}")
        json.dump(score_log, f, indent=4)

    logger.info("Combining complete")


_description_text = textwrap.dedent(
    """
    Combine multiple scoring files in PGS Catalog format (see 
    https://www.pgscatalog.org/downloads/ for details) to a 'long' table of columns 
    needed for variant matching and subsequent calculation.

    Custom scorefiles in PGS Catalog format can be combined with PGS Catalog scoring files, and 
    optionally liftover genomic coordinates to GRCh37 or GRCh38. The script can accept a mix of
    unharmonised and harmonised PGS Catalog data. By default all variants are output (including 
    positions with duplicated data [often caused by rsID/liftover collions across builds]) and 
    variants with missing positions. 
"""
)

_epilog_text = textwrap.dedent(
    """
    The long table is used to simplify intersecting variants in target genotyping datasets 
    and the scoring files with the match_variants program.
    """
)


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description=_description_text,
        epilog=_epilog_text,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-s",
        "--scorefiles",
        dest="scorefiles",
        nargs="+",
        help="<Required> Scorefile paths",
        required=True,
    )
    parser.add_argument(
        "--liftover",
        dest="liftover",
        help="<Optional> Convert scoring file variants to target genome build?",
        action="store_true",
    )
    parser.add_argument(
        "-t",
        "--target_build",
        dest="target_build",
        choices=["GRCh37", "GRCh38"],
        help="<Required> Build of target genome",
        required=True,
    )
    parser.add_argument(
        "-c",
        "--chain_dir",
        dest="chain_dir",
        help="Path to directory containing chain files",
        required="--liftover" in sys.argv,
    )
    parser.add_argument(
        "-m",
        "--min_lift",
        dest="min_lift",
        help="<Optional> If liftover, minimum proportion of variants lifted over",
        default=0.95,
        type=float,
    )
    parser.add_argument(
        "--drop_missing",
        dest="drop_missing",
        action="store_true",
        help="<Optional> Drop variants with missing information (chr/pos) and "
        "non-standard alleles (e.g. HLA=P/N) from the output file.",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        required=True,
        default="combined.txt",
        help="<Required> Output path to combined long scorefile "
        "[ will compress output if filename ends with .gz ]",
    )
    parser.add_argument(
        "-l",
        "--logfile",
        dest="logfile",
        default="log_combined.json",
        help="<Required> Name for the log file (score metadata) for combined scores."
        "[ will write to identical directory as combined scorefile]",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="<Optional> Extra logging information",
    )

    return parser.parse_args(args)


if __name__ == "__main__":
    run()
