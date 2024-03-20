import argparse
import logging
import pathlib
import textwrap
import operator
import functools

from ..lib.polygenicscore import PolygenicScore

logger = logging.getLogger(__name__)


def run_aggregate():
    logging.basicConfig(
        format="%(asctime)s %(name)s %(levelname)-8s %(message)s",
        level=logging.WARNING,
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    args = _parse_args()

    if args.verbose:
        logger.setLevel(logging.INFO)

    if not (outdir := pathlib.Path(args.outdir)).exists():
        raise FileNotFoundError(f"--outdir {outdir.name} doesn't exist")

    score_paths = [pathlib.Path(x) for x in args.scores]
    pgs = [PolygenicScore(path=x) for x in score_paths]
    # call __add__ a lot
    aggregated = functools.reduce(operator.add, pgs)
    aggregated.write(outdir=args.outdir, split=args.split)


def _description_text() -> str:
    return textwrap.dedent(
        """
    Aggregate plink .sscore files into a combined TSV table.

    This aggregation sums scores that were calculated from plink
    .scorefiles. Scorefiles may be split to calculate scores over different
    chromosomes or effect types. The PGS Catalog calculator automatically splits
    scorefiles where appropriate, and uses this script to combine them.

    Input .sscore files can be optionally compressed with zstd or gzip. 

    The aggregated output scores are compressed with gzip.
   """
    )


def _parse_args(args=None):
    parser = argparse.ArgumentParser(
        description=_description_text(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-s",
        "--scores",
        dest="scores",
        required=True,
        nargs="+",
        help="<Required> List of scorefile paths. Use a wildcard (*) to select multiple files.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        required=True,
        default="scores/",
        help="<Required> Output directory to store downloaded files",
    )
    parser.add_argument(
        "--split",
        dest="split",
        required=True,
        action=argparse.BooleanOptionalAction,
        help="Make one aggregated file per sampleset",
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
    run_aggregate()
