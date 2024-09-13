""" Public interface to the Polygenic Score Catalog core package """
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
    GenomeBuild,
    TargetVariants,
    TargetVariant,
    TargetType,
    RelabelArgs,
    relabel,
    relabel_write,
    effect_type_keyfunc,
    chrom_keyfunc,
    EffectType,
    models,
)

log_fmt = "%(name)s: %(asctime)s %(levelname)-8s %(message)s"
logging.basicConfig(format=log_fmt, datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)

__all__ = [
    "ScoringFiles",
    "ScoringFile",
    "Config",
    "GenomeBuild",
    "CatalogQuery",
    "ScoreQueryResult",
    "CatalogCategory",
    "TargetVariant",
    "TargetVariants",
    "TargetType",
    "NormalisedScoringFile",
    "RelabelArgs",
    "relabel",
    "relabel_write",
    "effect_type_keyfunc",
    "chrom_keyfunc",
    "EffectType",
    "models",
]

__version__ = importlib.metadata.version("pgscatalog.core")
