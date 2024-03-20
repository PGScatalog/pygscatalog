import argparse
import logging
import pathlib
import textwrap

from ..lib.polygenicscore import AdjustArguments, AggregatedPGS
from ..lib.principalcomponents import PrincipalComponents, PopulationType

logger = logging.getLogger(__name__)


def run_ancestry():
    logging.basicConfig(
        format="%(asctime)s %(name)s %(levelname)-8s %(message)s",
        level=logging.WARNING,
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    args = _parse_args()

    if args.verbose:
        logger.setLevel(logging.INFO)
        logger.info("Starting ancestry adjustment")
        logger.info("Verbose mode enabled")

    if not (outdir := pathlib.Path(args.outdir)).exists():
        raise FileNotFoundError(f"--outdir {outdir.name} doesn't exist")

    logger.info("Loading reference PCA data")
    ref_pc = PrincipalComponents(
        pcs_path=args.ref_pcs,
        dataset=args.d_ref,
        pop_type=PopulationType.REFERENCE,
        psam_path=args.psam,
        related_path=args.ref_related,
        npcs_popcomp=args.nPCs_popcomp,
        npcs_normalization=args.nPCs_normalization,
    )

    logger.info("Loading target PCA data")
    target_pcs = PrincipalComponents(
        pcs_path=args.target_pcs,
        dataset=args.d_target,
        pop_type=PopulationType.TARGET,
        npcs_popcomp=args.nPCs_popcomp,
        npcs_normalization=args.nPCs_normalization,
    )

    logger.info("Loading aggregated polygenic scores")
    # the PGS has already been aggregated and melted
    pgs = AggregatedPGS(path=args.scorefile, target_name=args.d_target)

    adjust_arguments = AdjustArguments(
        method_compare=args.method_compare,
        pThreshold=args.pThreshold,
        method_normalization=args.method_normalization,
    )
    logger.info(f"Set up parameters {adjust_arguments!r}")

    logger.info("Starting score adjustment")
    adjust_results = pgs.adjust(
        ref_pc=ref_pc, target_pc=target_pcs, adjust_arguments=adjust_arguments
    )

    logger.info(f"Writing results to {args.outdir}")
    adjust_results.write(directory=args.outdir)

    logger.info("Finished :]")


def _description_text() -> str:
    return textwrap.dedent(
        "Program to analyze ancestry outputs of the pgscatalog/pgsc_calc pipeline. Current inputs: "
        "\n  - PCA projections from reference and target datasets (*.pcs)"
        "\n  - calculated polygenic scores (e.g. aggregated_scores.txt.gz), "
        "\n  - information about related samples in the reference dataset (e.g. "
        "deg2_hg38.king.cutoff.out.id)."
    )


def _parse_args(args=None):
    parser = argparse.ArgumentParser(
        description=_description_text(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-d",
        "--dataset",
        dest="d_target",
        required=True,
        help="<Required> Label of the TARGET genomic dataset",
    )
    parser.add_argument(
        "-r",
        "--reference",
        dest="d_ref",
        required=True,
        help="<Required> Label of the REFERENCE genomic dataset",
    )
    parser.add_argument(
        "--ref_pcs",
        dest="ref_pcs",
        required=True,
        nargs="+",
        help="<Required> Principal components path (output from fraposa_pgsc)",
    )
    parser.add_argument(
        "--target_pcs",
        dest="target_pcs",
        required=True,
        nargs="+",
        help="<Required> Principal components path (output from fraposa_pgsc)",
    )
    parser.add_argument(
        "--psam",
        dest="psam",
        required=True,
        help="<Required> Reference sample information file path in plink2 psam format)",
    )
    parser.add_argument(
        "-x",
        "--reference_related",
        dest="ref_related",
        help="File of related sample IDs (excluded from training ancestry assignments)",
    )
    parser.add_argument(
        "-p",
        "--pop_label",
        dest="ref_label",
        default="SuperPop",
        help="Population labels in REFERENCE psam to use for assignment",
    )
    parser.add_argument(
        "-s",
        "--agg_scores",
        dest="scorefile",
        default="aggregated_scores.txt.gz",
        help="Aggregated scores in PGS Catalog format ([sampleset, IID] indexed)",
    )
    parser.add_argument(
        "-a",
        "--ancestry_method",
        dest="method_compare",
        choices=["RandomForest", "Mahalanobis"],
        default="RandomForest",
        help="Method used for population/ancestry assignment",
    )
    parser.add_argument(
        "--n_popcomp",
        dest="nPCs_popcomp",
        type=int,
        metavar="[1-20]",
        choices=range(1, 21),
        default=5,
        help="Number of PCs used for population comparison (default = 5)",
    )
    parser.add_argument(
        "-t",
        "--pval_threshold",
        dest="pThreshold",
        type=float,
        help="p-value threshold used to identify low-confidence ancestry similarities",
    )
    parser.add_argument(
        "-n",
        "--normalization_method",
        nargs="+",
        dest="method_normalization",
        choices=["empirical", "mean", "mean+var"],
        default=["empirical", "mean", "mean+var"],
        help="Method used for adjustment of PGS using genetic ancestry",
    )
    parser.add_argument(
        "--n_normalization",
        dest="nPCs_normalization",
        type=int,
        metavar="[1-20]",
        choices=range(1, 21),
        default=4,
        help="Number of PCs used for population NORMALIZATION (default = 4)",
    )
    parser.add_argument(
        "--outdir", dest="outdir", required=True, help="<Required> Output directory"
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
    run_ancestry()
