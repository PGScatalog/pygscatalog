import argparse
import concurrent.futures
import logging
import pathlib
import sys
import textwrap

from tqdm import tqdm

from pgscatalog.core.cli.normalise import write_normalised
from pgscatalog.core.lib.models import ScoreLog, ScoreLogs
from pgscatalog.core.lib import GenomeBuild, ScoringFile

logger = logging.getLogger(__name__)


def run():
    args = parse_args()
    if args.verbose:
        logging.getLogger("pgscatalog.core").setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    out_dir = pathlib.Path(args.outfile).resolve()

    if not out_dir.is_dir():
        raise NotADirectoryError

    paths = list(set(args.scorefiles))  # unique paths only
    scoring_files = sorted([ScoringFile(x) for x in paths], key=lambda x: x.pgs_id)
    target_build = GenomeBuild.from_string(args.target_build)

    for x in scoring_files:
        if x.genome_build is None and target_build is not None:
            raise ValueError(
                f"Can't combine {x.pgs_id} with missing build in "
                f"header when requesting {target_build=}"
            )

        if x.genome_build != target_build and not args.liftover:
            raise ValueError(
                f"Can't combine scoring file with genome build {x.genome_build!r} when {target_build=} without --liftover"
            )

    variant_log: list[ScoreLog] = []

    if args.liftover:
        chain_dir = pathlib.Path(args.chain_dir)
        if not chain_dir.exists():
            logger.critical(f"Chain directory is missing: {chain_dir}")
            raise FileNotFoundError

        liftover_kwargs = {
            "liftover": True,
            "chain_dir": chain_dir,
        }
    else:
        liftover_kwargs = {"liftover": False}

    n_finished = 0

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for scorefile in scoring_files:
            if scorefile.path.name.endswith(".gz"):
                gzip_output = True
            else:
                gzip_output = False

            out_path = out_dir / ("normalised_" + scorefile.path.name)

            futures.append(
                executor.submit(
                    write_normalised,
                    scorefile=scorefile,
                    target_build=target_build,
                    liftover_kwargs=liftover_kwargs,
                    drop_missing=args.drop_missing,
                    out_path=out_path,
                    gzip_output=gzip_output,
                    batch_size=args.batch_size,
                )
            )

        for future in tqdm(
            concurrent.futures.as_completed(futures), total=len(futures)
        ):
            log: ScoreLog = future.result()
            if log.compatible_effect_type:
                logger.info(
                    f"Finished scorefile with compatible effect type {log.pgs_id}"
                )
                n_finished += 1
            else:
                logger.info(
                    f"Couldn't process {log.pgs_id} because of incompatible effect type"
                )
            variant_log.append(log)

    if n_finished == 0 or not out_dir.exists():
        raise ValueError(
            "Couldn't process any scoring files. Did they all have non-additive weights?"
        )

    if n_finished != len(scoring_files):
        logger.warning(f"{len(scoring_files) - n_finished} scoring files were skipped")

    log_out_path = pathlib.Path(args.outfile) / args.logfile
    with open(log_out_path, "w") as f:
        logger.info(f"Writing log to {f.name}")
        f.write(ScoreLogs(variant_log).model_dump_json())

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
    parser.add_argument(
        "--batch_size",
        dest="batch_size",
        type=int,
        default=100_000,
        help="<Optional> Number of records to process in each batch",
    )
    parser.add_argument(
        "--threads",
        dest="threads",
        type=int,
        default=1,
        help="<Optional> Number of Python worker processes to use",
    )

    return parser.parse_args(args)


if __name__ == "__main__":
    run()
