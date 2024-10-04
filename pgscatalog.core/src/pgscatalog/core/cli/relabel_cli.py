import argparse
import logging
import pathlib

from pgscatalog.core.lib import relabel, relabel_write, RelabelArgs

logger = logging.getLogger(__name__)


def _parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Relabel the column values in one file based on a pair of columns in another",
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
        "-m",
        "--maps",
        help="mapping filenames",
        dest="map_files",
        nargs="+",
        required=True,
    )
    parser.add_argument(
        "-o", "--outdir", help="output directory", dest="outdir", required=True
    )
    parser.add_argument(
        "--col_from", help="column to change FROM", dest="col_from", required=True
    )
    parser.add_argument(
        "--col_to", help="column to change TO", dest="col_to", required=True
    )
    parser.add_argument(
        "--target_file", help="target file", dest="target_file", required=True
    )
    parser.add_argument(
        "--target_col",
        help="target column to revalue",
        dest="target_col",
        required=True,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="<Optional> Extra logging information",
    )
    parser.add_argument("--split", dest="split", action="store_true", required=False)
    parser.add_argument(
        "--combined", dest="combined", action="store_true", required=False
    )
    parser.add_argument("-cc", "--comment_char", dest="comment_char", default="##")
    args = parser.parse_args()

    if not (args.split or args.combined):
        parser.error("At least one of --combined or --split is required")

    return args


def run():
    args = _parse_args()

    if args.verbose:
        logging.getLogger("pgscatalog.corelib").setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    relabel_args = RelabelArgs(
        comment_char=args.comment_char,
        dataset=args.dataset,
        map_col_from=args.col_from,
        map_col_to=args.col_to,
        target_col=args.target_col,
    )
    logger.debug(f"Relabel arguments {relabel_args}")

    map_paths = [pathlib.Path(x) for x in args.map_files]
    in_path = pathlib.Path(args.target_file)

    for x in [*map_paths, in_path]:
        if not x.exists():
            raise FileNotFoundError(f"{x}")

    logger.debug("Relabelling variants")
    relabelled = relabel(
        in_path=in_path, map_paths=map_paths, relabel_args=relabel_args
    )
    logger.debug(f"Writing relabelled data to {args.outdir}")
    relabel_write(
        relabelled=relabelled,
        dataset=relabel_args.dataset,
        split_output=args.split,
        combined_output=args.combined,
        out_dir=args.outdir,
    )


if __name__ == "__main__":
    run()
