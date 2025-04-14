"""Convenience functions for pgscatalog-format"""

import csv
import logging
import pathlib
import shutil
import tempfile
from itertools import islice
from typing import Optional

import pydantic
from xopen import xopen

from pgscatalog.core import ScoringFile, GenomeBuild
from pgscatalog.core.lib import EffectTypeError
from pgscatalog.core.lib.models import ScoreLog, VariantLog


logger = logging.getLogger(__name__)


def batched(iterable, n):
    """
    Batch data into lists of length n. The last batch may be shorter.

    (stdlib itertools.batched was introduced in Python 3.12, so this function is helpful)
    """
    # batched('ABCDEFG', 3) --> ABC DEF G
    it = iter(iterable)
    while True:
        batch = tuple(islice(it, n))
        if not batch:
            return
        yield batch


def format_and_write(
    scorefile: ScoringFile,
    target_build: GenomeBuild,
    drop_missing: bool,
    out_path: pathlib.Path,
    liftover_kwargs: dict,
    batch_size: int = 100_000,
    gzip_output: bool = False,
) -> Optional[ScoreLog]:
    """Take a PGS Catalog scoring file and write a subset of fields to a consistent structure

    Returns a score log generated from the score header and variant statistics
    """
    logger.info(f"Started processing {scorefile.pgs_id}")
    variant_log: list[VariantLog] = []
    is_compatible: bool = True

    # suffix is important for xopen to automatically compress output
    suffix = ".txt.gz" if gzip_output else ".txt"

    temp_path = pathlib.Path(
        tempfile.NamedTemporaryFile(suffix=suffix, delete=False).name
    )

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
        # rename doesn't support cross-filesystem moves
        shutil.move(temp_path.resolve(), out_path.resolve())
        logger.info(f"Written normalised score to {out_path=}")
    finally:
        temp_path.unlink(missing_ok=True)

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
