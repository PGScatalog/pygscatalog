import importlib.metadata

from .lib import (
    VALID_CHROMOSOMES,
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
    TargetVariants,
)

__all__ = [
    "Scorefiles",
    "TargetGenome",
    "TargetVariants",
    "ScorePipeline",
    "GenomeFileType",
    "Pathish",
    "PathishList",
    "VALID_CHROMOSOMES",
    # legacy stuff
    "PolygenicScore",
    "PrincipalComponents",
    "PopulationType",
    "AggregatedPGS",
    "AdjustArguments",
    "AdjustResults",
]

__version__ = importlib.metadata.version("pgscatalog.calc")
