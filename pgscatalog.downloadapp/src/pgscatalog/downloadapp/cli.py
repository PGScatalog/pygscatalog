import argparse
import concurrent
import textwrap
from concurrent.futures import ThreadPoolExecutor

from pgscatalog.corelib import ScoringFiles, GenomeBuild
from pgscatalog.corelib import config


def run():
    args = parse_args()

    if args.pgsc_calc is not None:
        config.API_HEADER = {"user-agent": args.pgsc_calc}

    build = GenomeBuild.from_string(args.build)

    # unpack all accessions into a single flat list
    sfs = ScoringFiles([*args.pgs, *args.pgp, *args.efo], target_build=build)

    with ThreadPoolExecutor() as executor:
        futures = []
        for scorefile in sfs:
            futures.append(executor.submit(scorefile.download,
                                           overwrite=args.overwrite_existing_file,
                                           directory=args.outdir))

        for future in concurrent.futures.as_completed(futures):
            # nothing returned, but important to raise exceptions
            future.result()


description_text = textwrap.dedent(
    """\
Download a set of scoring files from the PGS Catalog using PGS
Scoring IDs, traits, or publication IDs.

The PGS Catalog API is queried to get a list of scoring file
URLs. Scoring files are downloaded via FTP to a specified
directory. PGS Catalog scoring files are staged with the name:

        {PGS_ID}.txt.gz

If a valid build is specified harmonized files are downloaded as:

    {PGS_ID}_hmPOS_{genome_build}.txt.gz

These harmonised scoring files contain genomic coordinates,
remapped from author-submitted information such as rsids.
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
        help="PGS Catalog ID(s) (e.g. PGS000001)"
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
        help="Download Harmonized Scores with Positions in Genome build: GRCh37 or "
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
        "--pgsc_calc",
        dest="pgsc_calc",
        help="<Optional> Provide information about downloading scoring files via "
             "pgsc_calc",
    )

    args = parser.parse_args(args)

    if args.pgs is None and args.efo is None and args.pgp is None:
        parser.error("Please provide at least one: --pgs, --efo, or --pgp")

    return args


if __name__ == "__main__":
    run()
