import logging
import importlib.metadata

from .lib import (
    Config,
    ScoreQueryResult,
    CatalogQuery,
    CatalogCategory,
    ScoringFiles,
    ScoringFile,
    NormalisedScoringFile,
    ScoreVariant,
    EffectType,
    GenomeBuild,
    TargetVariants,
    TargetVariant,
    TargetType,
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
    RelabelArgs,
    relabel,
    relabel_write,
    effect_type_keyfunc,
    chrom_keyfunc,
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
    "TargetType",
    "NormalisedScoringFile",
    "EffectType",
    "RelabelArgs",
    "relabel",
    "relabel_write",
    "effect_type_keyfunc",
    "chrom_keyfunc",
]

__version__ = importlib.metadata.version(__name__)
