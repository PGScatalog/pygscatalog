"""A declared sex, and what the calls say about it (ADR 0014)

A **Declared sex** is QC input only. It never derives ploidy, never fills a ploidy a
format could not express, and never reaches a dosage or a score - ADR 0007's rejection
of a pipeline-wide sex column as a *scoring* input stands in full. Everything here
compares a declaration against ploidy that was observed per call, and stops when the two
disagree.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pathlib

logger = logging.getLogger(__name__)

# plink and UK Biobank's coding. 0 is a legal value everywhere: it is how a user
# declares an unknown sex, and how they exempt a sample from the check (ADR 0014)
UNKNOWN, MALE, FEMALE = 0, 1, 2
VALID_SEX_CODES = frozenset({UNKNOWN, MALE, FEMALE})

SEX_MAP_COLUMNS = ("sample_id", "sex")


def read_sex_map(path: pathlib.Path) -> dict[str, int]:
    """Read a two-column sex map: sample_id, sex.

    Tab separated, with a header naming both columns. Codes are 0 unknown, 1 male,
    2 female. Any other value is refused naming the file, the line and the value - a
    mandatory declaration is worth nothing if the parser guesses at it (ADR 0014).

    ponytail: hand-rolled over open() rather than duckdb read_csv, which is how this
    repo reads its other tabular input. Two columns do not need a query engine, and
    duckdb's parse errors would have to be remapped to name the offending line anyway.
    """
    codes: dict[str, int] = {}

    with path.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        if tuple(header) != SEX_MAP_COLUMNS:
            raise ValueError(
                f"{path}: sex map header must be "
                f"{list(SEX_MAP_COLUMNS)}, found {header}"
            )

        # line 1 was the header, so data starts at 2 and the number in an error message
        # is the number the user's editor shows
        for line_number, line in enumerate(f, start=2):
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(SEX_MAP_COLUMNS):
                raise ValueError(
                    f"{path} line {line_number}: expected "
                    f"{len(SEX_MAP_COLUMNS)} tab separated columns, found {len(fields)}"
                )

            sample_id, code = fields[0], fields[1]

            try:
                sex = int(code)
            except ValueError:
                raise ValueError(
                    f"{path} line {line_number}: sex must be one of "
                    f"{sorted(VALID_SEX_CODES)}, found {code!r}"
                ) from None

            if sex not in VALID_SEX_CODES:
                raise ValueError(
                    f"{path} line {line_number}: sex must be one of "
                    f"{sorted(VALID_SEX_CODES)}, found {code!r}"
                )

            if sample_id in codes:
                raise ValueError(
                    f"{path} line {line_number}: sample {sample_id!r} is declared twice"
                )

            codes[sample_id] = sex

    if not codes:
        raise ValueError(f"{path}: sex map declares no samples")

    return codes


def check_sex_map_covers(
    sex_map: dict[str, int], samples: list[str], source: str
) -> None:
    """Every genotyped sample must carry a declaration.

    Coverage is asymmetric, and the asymmetry is forced by sample batching (ADR 0003).
    A genotyped sample missing from the map is a hard error, because a mandatory
    declaration that is silently incomplete is not a declaration, and those samples
    would otherwise pass the check unexamined with nothing recorded. A sample in the
    map that this shard does not hold is fine and expected: batches see a subset of the
    cohort, so one cohort-wide map is the normal way to run this, and rejecting extra
    rows would make the map a function of how the cohort was chopped up.
    """
    missing = [sample for sample in samples if sample not in sex_map]

    if missing:
        shown = ", ".join(missing[:10])
        suffix = "" if len(missing) <= 10 else f" (and {len(missing) - 10} more)"
        raise ValueError(
            f"{len(missing)} of {len(samples)} samples in {source} are missing from "
            f"the sex map: {shown}{suffix}. Every genotyped sample needs a "
            f"declaration. A map that misses all of them is usually keyed on a "
            f"different sample ID space."
        )
