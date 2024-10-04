"""This module contains functions for reading data from PGS Catalog files.
These functions aren't really meant to be imported outside corelib"""

import logging
from typing import Generator, Iterator

from xopen import xopen

from pgscatalog.core.lib.models import ScoreVariant

logger = logging.getLogger(__name__)


def read_rows_lazy(
    *,
    csv_reader: Iterator[list[str]],
    fields: list[str],
    name: str,
    wide: bool,
    row_nr: int,
) -> Generator[ScoreVariant, None, None]:
    """Read rows from an open scoring file and instantiate them as ScoreVariants"""
    for row in csv_reader:
        variant = dict(zip(fields, row))

        if wide:
            ew_col_idxs: list[int] = [
                i for i, x in enumerate(["effect_weight_" in x for x in fields]) if x
            ]
            for i, weight_name in zip(ew_col_idxs, [fields[i] for i in ew_col_idxs]):
                yield ScoreVariant(
                    **variant,  # type: ignore
                    **{  # type: ignore
                        "accession": weight_name,
                        "row_nr": row_nr,
                        "effect_weight": variant[weight_name],
                    },
                )
        else:
            yield ScoreVariant(**variant, **{"accession": name, "row_nr": row_nr})  # type: ignore

        row_nr += 1


def get_columns(path):
    """Grab column labels from a PGS Catalog scoring file. line_no is useful to skip the header"""
    with xopen(path, mode="rt") as f:
        for i, line in enumerate(f):
            if line.startswith("#"):
                continue
            line_no, cols = i, line.strip().split()
            if len(set(cols)) != len(cols):
                logger.critical(f"Duplicated column names: {cols}")
                raise ValueError

            return line_no, cols


def detect_wide(cols: list[str]) -> bool:
    """
    Check columns to see if multiple effect weights are present. Multiple effect weights must be present in the form:
    effect_weight_suffix1
    effect_weight_suffix2
    """
    if any(["effect_weight_" in x for x in cols]):
        logger.info("Wide scoring file detected with multiple effect weights")
        return True
    else:
        return False
