import logging

from .polygenicscore import (
    PolygenicScore,
    AggregatedPGS,
    AdjustResults,
    AdjustArguments,
)
from .principalcomponents import PopulationType, PrincipalComponents

log_fmt = "%(name)s: %(asctime)s %(levelname)-8s %(message)s"
logging.basicConfig(format=log_fmt, datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)

__all__ = [
    "PolygenicScore",
    "PrincipalComponents",
    "PopulationType",
    "AggregatedPGS",
    "AdjustArguments",
    "AdjustResults",
]
