from __future__ import annotations

import argparse
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any


def check_positive(value: Any) -> int:
    try:
        ivalue = int(value)
    except (TypeError, ValueError) as e:
        raise argparse.ArgumentTypeError(f"{value!r} is not a valid integer") from e

    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f"{ivalue} is not a valid positive integer")

    return ivalue


def absence_semantics(value: Any) -> str:
    """What it means for a variant to be absent from a target genome file.

    The two possibilities produce byte-identical files, so this can only ever be
    declared, never detected (ADR 0004). Only `unobserved` is implemented; `homref`
    is rejected loudly rather than silently ignored, because ignoring it would score
    a sparse-called callset as if every absent variant were simply unknown.
    """
    if value == "unobserved":
        return "unobserved"

    if value == "homref":
        raise argparse.ArgumentTypeError(
            "--absence_semantics 'homref' is not implemented yet. Absent variants are "
            "always treated as never assessed in this release."
        )

    raise argparse.ArgumentTypeError(
        f"{value!r} is not valid absence semantics, choose 'unobserved' or 'homref'"
    )


def zero_one_float(x: Any) -> float:
    try:
        floating = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"{x!r} not a floating-point literal"
        ) from None

    if floating <= 0.0 or floating > 1.0:
        raise argparse.ArgumentTypeError(f"{floating!r} not in range (0.0, 1.0]")
    return floating
