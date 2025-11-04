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
