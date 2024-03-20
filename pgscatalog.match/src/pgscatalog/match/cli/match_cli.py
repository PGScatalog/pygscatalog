import argparse
import atexit
import logging
import pathlib
import sys
import tempfile
import textwrap

import polars as pl

from ..lib import (
    VariantFrame,
    ScoringFileFrame,
    match_variants,
    MatchResult,
    MatchResults,
)

from ._config import Config
from ._write import write_matches

logger = logging.getLogger(__name__)


def _exit(cleanup):
    if cleanup:
        # should already be emptied by context managers
        Config.TMPDIR.rmdir()


def run_match():
    atexit.register(_exit, cleanup=Config.CLEANUP)
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

    Config.MATCHTMP = Config.OUTDIR / "matchtmp"
    Config.MATCHTMP.mkdir(exist_ok=False)
    Config.TMPDIR = Config.OUTDIR / "tmp"
    Config.TMPDIR.mkdir(exist_ok=False)
    Config.CHROM = args.chrom
    Config.MIN_OVERLAP = args.min_overlap

    if (n_target := len(args.target)) == 0:
        raise ValueError(f"{n_target=} must provide at least one target")

    Config.SPLIT = args.split
    Config.COMBINED = args.combined

    # parameters that control how the best match candidate is chosen
    # missing parameters will be set to defaults specified in matchlib
    Config.MATCH_PARAMS = {
        k: v for k in Config.MATCH_KWARGS if (v := getattr(args, k)) is not None
    }

    # start doing the work
    with ScoringFileFrame(
        path=Config.SCOREFILE,
        chrom=Config.CHROM,
        cleanup=Config.CLEANUP,
        tmpdir=Config.TMPDIR,
    ) as score_df, pl.StringCache():
        matchresults = []
        for target_f in args.target:
            ipc_path = get_match_candidates(
                target=target_f,
                score_df=score_df,
                dataset=Config.DATASET,
                tmpdir=Config.TMPDIR,
                chrom=Config.CHROM,
                cleanup=Config.CLEANUP,
            )
            # MatchResult is now backed by a stable-ish IPC file that will rarely get cleaned up
            matchresults.append(
                MatchResult.from_ipc(
                    matchresults_ipc_path=ipc_path, dataset=Config.DATASET
                )
            )

        if args.only_match:
            logger.warning("--only_match set, exiting with grace and aplomb")
            logger.warning(
                f"You'll need to combine match candidates using the Arrow IPC files in {Config.OUTDIR}"
            )
            sys.exit()
        else:
            matchresults = MatchResults(*matchresults)
            _ = write_matches(matchresults=matchresults, score_df=score_df)


def get_match_candidates(target, score_df, chrom, dataset, **kwargs):
    # don't clean this up, useful for debugging
    _, fout = tempfile.mkstemp(dir=Config.MATCHTMP)  # file names don't matter
    variants = VariantFrame(path=target, chrom=chrom, dataset=dataset, **kwargs)
    with variants as target_df:
        _ = match_variants(
            score_df=score_df, target_df=target_df, target=variants
        ).collect(outfile=fout)

    return fout


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description=_description_text(),
        epilog=_epilog_text(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-d",
        "--dataset",
        dest="dataset",
        required=True,
        help="<Required> Label for target genomic dataset",
    )
    parser.add_argument(
        "-s",
        "--scorefiles",
        dest="scorefile",
        required=True,
        help="<Required> Combined scorefile path (output of read_scorefiles.py)",
    )
    parser.add_argument(
        "-t",
        "--target",
        dest="target",
        required=True,
        nargs="+",
        help="<Required> A list of paths of target genomic variants (.bim format)",
    )
    parser.add_argument(
        "-c",
        "--chrom",
        dest="chrom",
        required=False,
        type=str,
        help="<Optional> Set which chromosome is in the target variant file to speed up matching ",
    )
    parser.add_argument(
        "--only_match",
        dest="only_match",
        action="store_true",
        help="<Optional> Only match, then write intermediate files, don't make scoring files",
    )
    parser.add_argument(
        "--min_overlap",
        dest="min_overlap",
        required=False,
        type=float,
        help="<Optional> Minimum proportion of variants to match before error",
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
        help="<Optional> Split scorefile per chromosome?",
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


def _description_text() -> str:
    return textwrap.dedent(
        """\
    Match variants from a combined scoring file against a set of
    target genomes from the same fileset, and output scoring files
    compatible with the plink2 --score function.

    A combined scoring file is the output of the combine_scorefiles
    script. It has the following structure:

        | chr_name | chr_position | ... | accession |
        | -------- | ------------ | --- | --------- |
        | 1        | 1            | ... | PGS000802 |

    The combined scoring file is in long format, with one row per
    variant for each scoring file (accession). This structure is
    different to the PGS Catalog standard, because the long format
    makes matching faster and simpler.

    Target genomes can be in plink1 bim format or plink2 pvar
    format. Variant IDs should be unique so that they can be specified
    in the scoring file as: variant_id|effect_allele|[effect_weight column(s)...] 

    Only one set of target genomes should be matched at a time. Don't
    try to match target genomes from different plink filesets. Matching 
    against a set of chromosomes from the same fileset is OK (see --split). 
   """
    )


def _epilog_text() -> str:
    return textwrap.dedent(
        """\
    match_variants will output at least one scoring file in a
    format compatible with the plink2 --score function. This
    output might be split across different files to ensure each
    variant ID, effect allele, and effect type appears only once
    in each file. Output files have the pattern:

        {dataset}_{chromosome}_{effect_type}_{n}.scorefile.

    If multiple chromosomes are combined into a single file (i.e. not
    --split), then {chromosome} is replaced with 'ALL'. Once the
    scorefiles are used to calculate a score with plink2, the .sscore
    files will need to be aggregated to calculate a single polygenic
    score for each dataset, sample, and accession (scoring file). The
    PGS Catalog Calculator does this automatically.
    """
    )


def _check_args(args):
    if args.chrom is not None and not args.only_match:
        # filtering the scoring file will break overlap assumptions and calculations
        # e.g.:
        #   what if one chromosome matches well but another chromosome matches poorly?
        #   what if the poorly matching chromosome only has 5 variants to match?
        #
        # pgsc_calc uses global overlap % to decide if a score fails matching
        # --only_match skips overlap calculations (done in combine_matches instead)
        logger.critical("--chrom requires --only_match")
        sys.exit(1)
    if args.only_match and args.min_overlap is not None:
        # can't calculate min_overlap properly if just checking matches
        logger.critical("Invalid arguments: --only_match and --min_overlap (pick one!)")
        sys.exit(1)
    if not args.only_match and args.min_overlap is None:
        # need to calculate min_overlap before making scoring files
        logger.critical("Invalid arguments: set --min_overlap or --only_match")
        sys.exit(1)
    if args.split and args.only_match:
        # not writing scoring files, so split output doesn't make sense
        logger.critical("Invalid arguments: --only_match and --split (pick one!)")
        sys.exit(1)
    if any(
        [
            x in sys.argv
            for x in [
                "--keep_first_match",
                "--ignore_strand_flips",
                "--keep_multiallelic",
                "--keep_ambiguous",
            ]
        ]
    ):
        logger.warning(
            "Invalid arguments: --only_match and --keep_first_match, --ignore_strand_flips,"
            "keep_multiallelic, or keep_ambiguous"
        )
        logger.warning("Pass these arguments to combine_matches instead")

    if args.combined is False and args.split is False:
        logger.warning("No output format specified, writing to combined scoring file")
        args.combined = True

    return args


if __name__ == "__main__":
    run_match()
