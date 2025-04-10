"""This module contains functions useful for normalising scoring files into a
standardised format, combining them, and calculating some statistics. Only really
useful for the CLI, not good for importing elsewhere."""

import csv
import functools
import gzip
import logging
import os


logger = logging.getLogger(__name__)


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
            self.open_function = open

    def write(self, batch):
        mode = "at" if os.path.exists(self.filename) else "wt"
        with self.open_function(self.filename, mode) as f:
            writer = csv.DictWriter(
                f,
                delimiter="\t",
                lineterminator="\n",
                fieldnames=self.fieldnames,
                extrasaction="ignore",
            )
            match mode:
                case "wt":
                    logger.info("Writing variants to new file")
                    writer.writeheader()
                    writer.writerows(batch)
                case "at":
                    logger.info("Appending variants to existing file")
                    writer.writerows(batch)
                case _:
                    raise ValueError(f"Invalid {mode=}")
