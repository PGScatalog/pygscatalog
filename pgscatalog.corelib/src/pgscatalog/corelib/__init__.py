import logging

from ._config import Config
from .catalogapi import ScoreQueryResult, CatalogQuery, CatalogCategory
from .scorefiles import ScoringFiles, ScoringFile, NormalisedScoringFile
from .scorevariant import ScoreVariant, EffectType
from .genomebuild import GenomeBuild
from .targetvariants import TargetVariants, TargetVariant
from .pgsexceptions import (
    BasePGSException,
    MatchError,
    DuplicateMatchError,
    MatchRateError,
    ZeroMatchesError,
    MatchValueError,
    CombineError,
    BuildError,
    ScoreFormatError,
    CatalogError,
    ScoreDownloadError,
    ScoreChecksumError,
    QueryError,
    InvalidAccessionError,
    SamplesheetError,
    GenomesNotFound,
    SamplesheetFormatError,
)

log_fmt = "%(name)s: %(asctime)s %(levelname)-8s %(message)s"
logging.basicConfig(format=log_fmt, datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)

__all__ = [
    "BasePGSException",
    "MatchError",
    "DuplicateMatchError",
    "MatchRateError",
    "ZeroMatchesError",
    "MatchValueError",
    "CombineError",
    "BuildError",
    "ScoreFormatError",
    "CatalogError",
    "ScoreDownloadError",
    "ScoreChecksumError",
    "QueryError",
    "InvalidAccessionError",
    "SamplesheetError",
    "GenomesNotFound",
    "SamplesheetFormatError",
    "ScoringFiles",
    "ScoringFile",
    "ScoreVariant",
    "Config",
    "GenomeBuild",
    "CatalogQuery",
    "ScoreQueryResult",
    "CatalogCategory",
    "TargetVariant",
    "TargetVariants",
    "NormalisedScoringFile",
    "EffectType",
]
