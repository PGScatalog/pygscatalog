import argparse
import concurrent
import logging
import pathlib
import textwrap
from concurrent.futures import ThreadPoolExecutor

from tqdm import tqdm

from ..lib import ScoringFiles, GenomeBuild, Config

logger = logging.getLogger(__name__)


def run():
    args = parse_args()

    if args.verbose:
        logging.getLogger("pgscatalog.corelib").setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    if not (outdir := pathlib.Path(args.outdir)).exists():
        raise FileNotFoundError(f"--outdir {outdir.name} doesn't exist")

    if args.user_agent is not None:
        logger.info(f"Setting user agent to {args.user_agent}")
        Config.API_HEADER = {"user-agent": args.user_agent}

    build = GenomeBuild.from_string(args.build)
    if build is None:
        logger.warning(
            "--build is missing, downloading scoring files in author-reported build"
        )
    else:
        logger.info(f"Downloading scoring files that have been harmonised to {build=}")

    # unpack all accessions into a single flat list
    sfs = ScoringFiles(
        [*args.pgs, *args.pgp, *args.efo],
        target_build=build,
        include_children=args.efo_include_children,
    )

    with ThreadPoolExecutor() as executor:
        futures = []
        for scorefile in sfs:
            logger.info(f"Submitting {scorefile!r} download")
            futures.append(
                executor.submit(
                    scorefile.download,
                    overwrite=args.overwrite_existing_file,
                    directory=args.outdir,
                )
            )

        for future in tqdm(
            concurrent.futures.as_completed(futures), total=len(futures)
        ):
            # nothing returned, but important to grab result to raise exceptions
            future.result()
            logger.info("Download complete")

    logger.info("All downloads finished")


description_text = textwrap.dedent(
    """\
Download a set of scoring files from the PGS Catalog using PGS Scoring
IDs, traits, or publication accessions.

The PGS Catalog API is queried to get a list of scoring file URLs.
Scoring files are downloaded asynchronously via HTTPS to a specified
directory. Downloaded files are automatically validated against an md5
checksum.

PGS Catalog scoring files are staged with the name:

    {PGS_ID}.txt.gz

If a valid build is specified harmonized files are downloaded as:

    {PGS_ID}_hmPOS_{genome_build}.txt.gz

These harmonised scoring files contain genomic coordinates, remapped
from author-submitted information such as rsIDs.
"""
)


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description=description_text,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--pgs",
        nargs="+",
        dest="pgs",
        default=[],
        help="PGS Catalog ID(s) (e.g. PGS000001)",
    )
    parser.add_argument(
        "-t",
        "--efo",
        dest="efo",
        default=[],
        nargs="+",
        help="Traits described by an EFO term(s) (e.g. EFO_0004611)",
    )
    parser.add_argument(
        "-e",
        "--efo_direct",
        dest="efo_include_children",
        action="store_false",
        help="<Optional> Return only PGS tagged with exact EFO term "
        "(e.g. no PGS for child/descendant terms in the ontology)",
    )
    parser.add_argument(
        "-p",
        "--pgp",
        dest="pgp",
        default=[],
        help="PGP publication ID(s) (e.g. PGP000007)",
        nargs="+",
    )
    parser.add_argument(
        "-b",
        "--build",
        dest="build",
        choices=["GRCh37", "GRCh38"],
        help="Download harmonized scores with positions in genome build: GRCh37 or "
        "GRCh38",
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
        "-w",
        "--overwrite",
        dest="overwrite_existing_file",
        action="store_true",
        help="<Optional> Overwrite existing Scoring File if a new version is "
        "available for download on the FTP",
    )
    parser.add_argument(
        "-c",
        "--user_agent",
        dest="user_agent",
        help="<Optional> Provide custom user agent when querying PGS Catalog API",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="<Optional> Extra logging information",
    )

    args = parser.parse_args(args)

    if args.pgs is None and args.efo is None and args.pgp is None:
        parser.error("Please provide at least one: --pgs, --efo, or --pgp")

    return args


if __name__ == "__main__":
    run()
