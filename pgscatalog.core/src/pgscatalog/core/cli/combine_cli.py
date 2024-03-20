import argparse
import concurrent.futures
import json
import logging
import pathlib
import sys
import textwrap

from tqdm import tqdm
from ..lib import GenomeBuild, ScoringFile

from ._combine import normalise, get_variant_log, TextFileWriter

logger = logging.getLogger(__name__)


def run():
    args = parse_args()

    if args.verbose:
        logging.getLogger("pgscatalog.corelib").setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    out_path = pathlib.Path(args.outfile)

    if out_path.exists():
        logger.critical(f"Output file already exists: {args.outfile}")
        raise FileExistsError

    match x := out_path.name:
        case _ if x.endswith("gz"):
            logger.debug("Compressing output with gzip")
            compress_output = True
        case _:
            logger.debug("Not compressing output")
            compress_output = False

    paths = list(set(args.scorefiles))  # unique paths only
    scoring_files = [ScoringFile(x) for x in paths]
    target_build = GenomeBuild.from_string(args.target_build)

    for x in scoring_files:
        if x.genome_build != target_build and not args.liftover:
            raise ValueError(
                f"Can't combine scoring file with genome build {x.genome_build!r} when {target_build=} without --liftover"
            )

    variant_log = []

    if args.liftover:
        chain_dir = pathlib.Path(args.chain_dir)
        if not chain_dir.exists():
            logger.critical(f"Chain directory is missing: {chain_dir}")
            raise FileNotFoundError

        liftover_kwargs = {
            "liftover": True,
            "chain_dir": chain_dir,
            "target_build": target_build,
        }
    else:
        liftover_kwargs = {"liftover": False}

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []
        for scorefile in scoring_files:
            logger.info(f"Submitting {scorefile!r} for execution")
            futures.append(
                executor.submit(
                    normalise,
                    scorefile,
                    drop_missing=args.drop_missing,
                    **liftover_kwargs,
                )
            )
        # TODO: do this I/O asynchronously
        for future in tqdm(
            concurrent.futures.as_completed(futures), total=len(futures)
        ):
            writer = TextFileWriter(compress=compress_output, filename=out_path)
            normalised_score = future.result()
            writer.write(normalised_score)
            variant_log.append(get_variant_log(normalised_score))

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
