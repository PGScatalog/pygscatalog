"""Public interface to the Polygenic Score Catalog core package"""

import importlib.metadata
import logging

from pgscatalog.core.lib import (
    CatalogCategory,
    CatalogQuery,
    Config,
    EffectType,
    GenomeBuild,
    RelabelArgs,
    ScoreQueryResult,
    ScoringFile,
    ScoringFiles,
    chrom_keyfunc,
    effect_type_keyfunc,
    models,
    relabel,
    relabel_write,
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
    "RelabelArgs",
    "relabel",
    "relabel_write",
    "effect_type_keyfunc",
    "chrom_keyfunc",
    "EffectType",
    "models",
]

__version__ = importlib.metadata.version("pgscatalog.core")
