import argparse
import csv
import logging
import pathlib
import textwrap
import typing
from collections import deque
from typing import Optional
from xopen import xopen

from ..lib import PolygenicScore
from pgscatalog.core import chrom_keyfunc

logger = logging.getLogger(__name__)


def get_id_from_scorefile(path: pathlib.Path) -> typing.Generator[str, None, None]:
    """Get IDs from a scoring file, where ID is the first column"""
    with xopen(path) as f:
        x = csv.reader(f, delimiter="\t")
        for i, field in enumerate(x):
            if i == 0:
                continue  # skip header
            else:
                yield field[0].strip()


def verify_variants(score: pathlib.Path) -> None:
    """Check that the variants used to calculate the score (.vars file) are identical to the variants in the scoring file (.scorefile.gz)
    Raises a ValueError if the variants don't overlap perfectly"""
    base_path: str = score.stem.split(".")[0]
    var_path: pathlib.Path = score.parent / (base_path + ".sscore.vars")
    scorefile_path: pathlib.Path = score.parent / (base_path + ".scorefile.gz")

    if not var_path.exists():
        raise FileNotFoundError(f"{var_path} doesn't exist")

    if not scorefile_path.exists():
        raise FileNotFoundError(f"{scorefile_path} doesn't exist")

    with open(var_path) as f:
        # just a file where each line contains a variant ID
        scored_variants = frozenset(x.strip() for x in f)
        logger.info(f"Read {len(scored_variants)} from {var_path}")

    scorefile_variants = frozenset(get_id_from_scorefile(scorefile_path))
    logger.info(f"Read {len(scorefile_variants)} from {scorefile_path}")

    if diff := scored_variants.symmetric_difference(scorefile_variants):
        raise ValueError(f"Missing variants {diff}")
    else:
        logger.info(f"{score} was calculated from good variant sets :)")


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

    if args.verify_variants:
        logger.info("Checking variant overlap")
        [verify_variants(x) for x in score_paths]

    # dfs are only read into memory after accessing them explicitly e.g. pgs[0].df
    pgs = deque(PolygenicScore(path=x) for x in score_paths)

    observed_columns = set()
    aggregated: Optional[PolygenicScore] = None

    # first, use PolygenicScore's __add__ method, which implements df.add(fill_value=0)
    while pgs:
        # popleft ensures that dfs are removed from memory after each aggregation
        score: PolygenicScore = pgs.popleft()
        if aggregated is None:
            logger.info(f"Initialising aggregation with {score}")
            aggregated: PolygenicScore = score
        else:
            logger.info(f"Adding {score}")
            aggregated += score
        observed_columns.update(set(score.df.columns))

    # check to make sure that every column we saw in the dataframes is in the output
    if (dfcols := set(aggregated.df.columns)) != observed_columns:
        raise ValueError(
            f"Missing columns in aggregated file!. "
            f"Observed: {observed_columns}. "
            f"In aggregated: {dfcols}"
        )
    else:
        logger.info("Aggregated columns match observed columns")

    # next, melt the plink2 scoring files from wide (many columns) format to long format
    aggregated.melt()

    # recalculate PGS average using aggregated SUM and DENOM
    aggregated.average()

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
    parser.add_argument(
        "--verify_variants",
        dest="verify_variants",
        action="store_true",
        help="<Optional> Verify variants from scoring file match scored variants perfectly."
        " Note: requires .scorefile.gz files used by plink to create the calculated score, and a .sscore.vars file output by plink after scoring"
        " It's assumed that these files have a common file name prefix (this is default plink behaviour) e.g. cineca_22_additive_0.sscore.zst needs cineca_22_additive_0.scorefile.gz & cineca_22_additive_0.sscore.vars",
    )
    return parser.parse_args(args)


if __name__ == "__main__":
    run_aggregate()
