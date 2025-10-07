import importlib.metadata

from .lib import (
    AdjustArguments,
    AdjustResults,
    AggregatedPGS,
    GenomeFileType,
    Pathish,
    PathishList,
    # legacy stuff
    PolygenicScore,
    PopulationType,
    PrincipalComponents,
    Scorefiles,
    ScorePipeline,
    TargetGenome,
    TargetVariant,
    TargetVariants,
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
