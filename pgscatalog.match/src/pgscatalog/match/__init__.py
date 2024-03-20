import logging

from .lib.variantframe import VariantFrame
from .lib.scoringfileframe import ScoringFileFrame, match_variants
from .lib.matchresult import MatchResult, MatchResults
from .lib.plinkscorefiles import PlinkScoreFiles


log_fmt = "%(name)s: %(asctime)s %(levelname)-8s %(message)s"
logging.basicConfig(format=log_fmt, datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)

__all__ = [
    "VariantFrame",
    "ScoringFileFrame",
    "MatchResult",
    "MatchResults",
    "PlinkScoreFiles",
    "match_variants",
]
