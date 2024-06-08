import argparse
import logging
import pathlib
import textwrap
from collections import deque

from ..lib import PolygenicScore
from pgscatalog.core import chrom_keyfunc

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
        logging.getLogger("pgscatalog.core").setLevel(logging.INFO)
        logging.getLogger("pgscatalog.calc").setLevel(logging.INFO)

    if not (outdir := pathlib.Path(args.outdir)).exists():
        raise FileNotFoundError(f"--outdir {outdir.name} doesn't exist")

    score_paths = sorted([pathlib.Path(x) for x in args.scores], key=chrom_keyfunc())
    # dfs are only read into memory after accessing them explicitly e.g. pgs[0].df
    pgs = deque(PolygenicScore(path=x) for x in score_paths)

    observed_columns = set()
    aggregated = None

    while pgs:
        x = pgs.popleft()  # to remove dfs from memory during iteration
        if aggregated is None:
            logger.info(f"Initialising aggregation with {x}")
            aggregated = x
        else:
            logger.info(f"Adding {x}")
            aggregated += x
        observed_columns.update(set(x.df.columns))

    if (dfcols := set(aggregated.df.columns)) != observed_columns:
        raise ValueError(
            f"Missing columns in aggregated file!. "
            f"Observed: {observed_columns}. "
            f"In aggregated: {dfcols}"
        )
    else:
        logger.info("Aggregated columns match observed columns")

    logger.info("Aggregation finished! Writing to a file")
    aggregated.write(outdir=args.outdir, split=args.split)
    logger.info("all done. bye :)")


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
