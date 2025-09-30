import importlib.metadata

from .lib import (
    Scorefiles,
    TargetGenome,
    TargetVariant,
    TargetVariants,
    ScorePipeline,
    GenomeFileType,
    Pathish,
    PathishList,
    # legacy stuff
    PolygenicScore,
    AggregatedPGS,
    AdjustResults,
    AdjustArguments,
    PopulationType,
    PrincipalComponents,
)

__all__ = [
    "Scorefiles",
    "TargetGenome",
    "TargetVariant",
    "TargetVariants",
    "ScorePipeline",
    "GenomeFileType",
    "Pathish",
    "PathishList",
    # legacy stuff
    "PolygenicScore",
    "PrincipalComponents",
    "PopulationType",
    "AggregatedPGS",
    "AdjustArguments",
    "AdjustResults",
]

__version__ = importlib.metadata.version("pgscatalog.calc")
