from .genomefiletypes import GenomeFileType

# legacy stuff
from .legacy.polygenicscore import (
    AdjustArguments,
    AdjustResults,
    AggregatedPGS,
    PolygenicScore,
)
from .legacy.principalcomponents import PopulationType, PrincipalComponents
from .scorefile import Scorefiles
from .scorepipeline import ScorePipeline
from .targetgenome import TargetGenome
from .targetvariants import TargetVariants
from .types import Pathish, PathishList

__all__ = [
    "Scorefiles",
    "TargetGenome",
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
    "PrincipalComponents",
]
