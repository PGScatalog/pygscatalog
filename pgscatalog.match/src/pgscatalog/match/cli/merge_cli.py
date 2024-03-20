import argparse
import logging
import pathlib

from ._config import Config
from ._write import write_matches

from ..lib import ScoringFileFrame, MatchResult, MatchResults

import polars as pl

logger = logging.getLogger(__name__)


def run_merge():
    args = parse_args()

    if args.verbose:
        logging.getLogger("pgscatalog.core").setLevel(logging.DEBUG)
        logging.getLogger("pgscatalog.match").setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    Config.DATASET = args.dataset
    Config.CLEANUP = False
    Config.OUTDIR = pathlib.Path(args.outdir)
    Config.SCOREFILE = pathlib.Path(args.scorefile)

    if not Config.OUTDIR.exists():
        raise FileNotFoundError(f"{Config.OUTDIR} does not exist")

    Config.TMPDIR = Config.OUTDIR / "tmp"
    Config.TMPDIR.mkdir(exist_ok=False)

    # parameters that control how the best match candidate is chosen
    # missing parameters will be set to defaults specified in matchlib
    Config.MATCH_PARAMS = {
        k: v for k in Config.MATCH_KWARGS if (v := getattr(args, k)) is not None
    }
    Config.MIN_OVERLAP = args.min_overlap
    Config.SPLIT = args.split
    Config.COMBINED = args.combined

    with ScoringFileFrame(
        path=Config.SCOREFILE,
        chrom=None,  # when merging, scoring files can't be filtered
        cleanup=Config.CLEANUP,
        tmpdir=Config.TMPDIR,
    ) as score_df, pl.StringCache():
        matchresults = MatchResults(
            *(MatchResult.from_ipc(x, dataset=Config.DATASET) for x in args.matches)
        )
        matchdf = write_matches(matchresults=matchresults, score_df=score_df)
        _check_duplicate_vars(matchdf)


def _check_duplicate_vars(matches):
    max_occurrence = (
        matches.filter(pl.col("match_status") == "matched")
        .group_by(["accession", "ID"])
        .len()
        .select("len")
        .max()
        .collect()
        .item(0, 0)
    )

    match n := max_occurrence:
        case None:
            logger.critical("No variant matches found")
            logger.critical(
                "Did you set the correct genome build? Did you impute your genomes?"
            )
            raise ValueError
        case _ if n > 1:
            logger.critical("Duplicate IDs in final matches")
            logger.critical(
                "Please double check your genomes for duplicates and try again"
            )
            raise ValueError
        case _:
            logger.info("Scoring files are valid (no duplicate variants found)")


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--dataset",
        dest="dataset",
        required=True,
        help="<Required> Label for target genomic dataset",
    )
    parser.add_argument(
        "-s",
        "--scorefile",
        dest="scorefile",
        required=True,
        help="<Required> Path to scorefile",
    )
    parser.add_argument(
        "-m",
        "--matches",
        dest="matches",
        required=True,
        nargs="+",
        help="<Required> List of match files",
    )
    parser.add_argument(
        "--min_overlap",
        dest="min_overlap",
        required=True,
        type=float,
        default=0.75,
        help="<Required> Minimum proportion of variants to match before error",
    )
    parser.add_argument(
        "-IDs",
        "--filter_IDs",
        dest="filter_IDs",
        help="<Optional> Path to file containing list of variant IDs that can be included in the final scorefile."
        "[useful for limiting scoring files to variants present in multiple datasets]",
    )
    parser.add_argument(
        "--outdir", dest="outdir", required=True, help="<Required> Output directory"
    )
    parser.add_argument(
        "--split",
        dest="split",
        default=False,
        action="store_true",
        help="<Optional> Write scorefiles split per chromosome?",
    )
    parser.add_argument(
        "--combined",
        dest="combined",
        default=False,
        action="store_true",
        help="<Optional> Write scorefiles in combined format?",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="<Optional> Extra logging information",
    )
    # variant matching arguments -------------------------------------------------------
    parser.add_argument(
        "--keep_ambiguous",
        dest="remove_ambiguous",
        action="store_false",
        help="""<Optional> Flag to force the program to keep variants with
                        ambiguous alleles, (e.g. A/T and G/C SNPs), which are normally
                        excluded (default: false). In this case the program proceeds
                        assuming that the genotype data is on the same strand as the
                        GWAS whose summary statistics were used to construct the score.
                                """,
    )
    parser.add_argument(
        "--keep_multiallelic",
        dest="remove_multiallelic",
        action="store_false",
        help="<Optional> Flag to allow matching to multiallelic variants (default: false).",
    )
    parser.add_argument(
        "--ignore_strand_flips",
        dest="skip_flip",
        action="store_true",
        help="""<Optional> Flag to not consider matched variants that may be reported 
                        on the opposite strand.  Default behaviour is to flip/complement unmatched variants and check if
                        they match.""",
    )
    parser.add_argument(
        "--keep_first_match",
        dest="keep_first_match",
        action="store_true",
        help="""<Optional> If multiple match candidates for a variant exist that can't be prioritised,
                         keep the first match candidate (default: drop all candidates)""",
    )

    return _check_args(parser.parse_args(args))


def _check_args(args):
    if args.combined is False and args.split is False:
        logger.warning("No output format specified, writing to combined scoring file")
        args.combined = True

    return args


if __name__ == "__main__":
    run_merge()
