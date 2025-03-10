import argparse
import concurrent.futures
import csv
import logging
import pathlib
import sys
import tempfile
import textwrap

import pydantic
from tqdm import tqdm
from xopen import xopen

from itertools import islice
from pgscatalog.core.lib.models import ScoreLog, ScoreLogs, VariantLog
from pgscatalog.core.lib import GenomeBuild, ScoringFile, EffectTypeError


logger = logging.getLogger(__name__)


def batched(iterable, n):
    """
    Batch data into lists of length n. The last batch may be shorter.

    itertools.batched was introduced in Python 3.12
    """
    # batched('ABCDEFG', 3) --> ABC DEF G
    it = iter(iterable)
    while True:
        batch = tuple(islice(it, n))
        if not batch:
            return
        yield batch


def _combine(
    scorefile: ScoringFile,
    target_build: GenomeBuild,
    drop_missing: bool,
    out_path: pathlib.Path,
    liftover_kwargs: dict,
    batch_size: int = 100_000,
) -> ScoreLog:
    """This function normalises a single scoring file to a consistent structure and returns a summary log generated from the score header and variant statistics"""
    logger.info(f"Started processing {scorefile.pgs_id}")
    variant_log: list[VariantLog] = []
    is_compatible: bool = True

    # suffix is important for xopen to compress output
    with tempfile.NamedTemporaryFile(suffix=".txt.gz") as temp:
        temp_path = pathlib.Path(temp.name)
        normalised_score = scorefile.normalise(
            drop_missing=drop_missing,
            **liftover_kwargs,
            target_build=target_build,
        )
        score_batches = batched(normalised_score, n=batch_size)

        try:
            with xopen(temp_path, mode="w") as out_csv:
                logger.info("Opening file for writing")
                output_field_order = [
                    "chr_name",
                    "chr_position",
                    "effect_allele",
                    "other_allele",
                    "effect_weight",
                    "effect_type",
                    "is_duplicated",
                    "accession",
                    "row_nr",
                ]
                writer = csv.DictWriter(
                    out_csv,
                    fieldnames=output_field_order,
                    delimiter="\t",
                    extrasaction="ignore",
                )
                writer.writeheader()
                for batch in score_batches:
                    variant_batch = [x.model_dump() for x in batch]
                    variant_log.extend([VariantLog(**x) for x in variant_batch])
                    logger.info(f"Writing variant batch (size: {batch_size}) to file")
                    writer.writerows(variant_batch)
        except EffectTypeError:
            logger.warning(
                f"Unsupported non-additive effect types in {scorefile=}, skipping"
            )
            is_compatible = False
        except pydantic.ValidationError:
            logger.critical(
                f"{scorefile.pgs_id} contains invalid data, stopping and exploding"
            )
            raise
        else:
            logger.info(f"Finished processing {scorefile.pgs_id}")
            temp_path.rename(out_path)
            logger.info(f"Written normalised score to {out_path=}")

    log: ScoreLog = ScoreLog(
        header=scorefile.header,
        variant_logs=variant_log,
        compatible_effect_type=is_compatible,
    )

    if log.variants_are_missing:
        logger.warning(
            f"{log.variant_count_difference} fewer variants in output compared to original file"
        )

    return log


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

    # TODO: max_workers=1: ready to parallelise and write to separate files
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        futures = []
        for scorefile in scoring_files:
            out_path = out_dir / ("normalised_" + scorefile.path.name)
            futures.append(
                executor.submit(
                    _combine,
                    scorefile=scorefile,
                    target_build=target_build,
                    liftover_kwargs=liftover_kwargs,
                    drop_missing=args.drop_missing,
                    out_path=out_path,
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

    return parser.parse_args(args)


if __name__ == "__main__":
    run()
