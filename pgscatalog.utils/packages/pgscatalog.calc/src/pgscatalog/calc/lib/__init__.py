from .genomefiletypes import GenomeFileType
from .scorefile import Scorefiles
from .scorepipeline import ScorePipeline
from .targetgenome import TargetGenome
from .targetvariant import TargetVariant, TargetVariants
from .types import Pathish, PathishList

# legacy stuff
from .legacy.polygenicscore import (
    PolygenicScore,
    AggregatedPGS,
    AdjustResults,
    AdjustArguments,
)
from .legacy.principalcomponents import PopulationType, PrincipalComponents

__all__ = [
    "Scorefiles",
    "TargetGenome",
    "TargetVariant",
    "TargetVariants",
    "ScorePipeline",
    "GenomeFileType",
    "Pathish",
    "PathishList",
    # legacy stuff exported
    "PolygenicScore",
    "PopulationType",
    "AggregatedPGS",
    "AdjustResults",
    "AdjustArguments",
    "PrincipalComponents"
]
