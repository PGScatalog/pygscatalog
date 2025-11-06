from .cache.genomefiletypes import GenomeFileType
from .cache.targetgenome import TargetGenome
from .cache.targetvariants import TargetVariants

# legacy stuff
from .legacy.polygenicscore import (
    AdjustArguments,
    AdjustResults,
    AggregatedPGS,
    PolygenicScore,
)
from .legacy.principalcomponents import PopulationType, PrincipalComponents
from .score.scorepipeline import ScorePipeline
from .scorefile import Scorefiles
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
