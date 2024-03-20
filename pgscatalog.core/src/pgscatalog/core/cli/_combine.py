""" This module contains functions useful for normalising scoring files into a
standardised format, combining them, and calculating some statistics. Only really
useful for the CLI, not good for importing elsewhere."""
import collections
import csv
import functools
import gzip
import logging
import os

from ..lib.scorevariant import ScoreVariant

logger = logging.getLogger(__name__)


def normalise(
    scoring_file, liftover=False, drop_missing=False, chain_dir=None, target_build=None
):
    """Normalise all rows in a scoring file"""
    return list(
        scoring_file.normalise(
            liftover=liftover,
            drop_missing=drop_missing,
            chain_dir=chain_dir,
            target_build=target_build,
        )
    )


def get_variant_log(batch):
    # these statistics can only be generated while iterating through variants
    n_variants = collections.Counter("n_variants" for item in batch)
    hm_source = collections.Counter(getattr(item, "hm_source") for item in batch)
    return n_variants + hm_source


class DataWriter:
    def __init__(self, filename):
        self.filename = filename
        self.fieldnames = [
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
        logger.info(f"Output filename: {filename}")

    def write(self, batch):
        pass


class TextFileWriter(DataWriter):
    def __init__(self, compress, filename):
        super().__init__(filename)
        self.compress = compress

        if self.compress:
            logger.info("Writing with gzip")
            self.open_function = functools.partial(gzip.open, compresslevel=6)
        else:
            logger.info("Writing text file")
            self.open_function = open

    def write(self, batch):
        mode = "at" if os.path.exists(self.filename) else "wt"
        with self.open_function(self.filename, mode) as f:
            writer = csv.writer(
                f,
                delimiter="\t",
                lineterminator="\n",
            )
            if mode == "wt":
                writer.writerow(ScoreVariant.output_fields)

            writer.writerows(batch)
