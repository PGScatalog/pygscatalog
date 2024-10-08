import logging

from pgscatalog.match.lib.variantframe import VariantFrame
from pgscatalog.match.lib.scoringfileframe import ScoringFileFrame, match_variants
from pgscatalog.match.lib.matchresult import MatchResult, MatchResults
from pgscatalog.match.lib.plinkscorefiles import PlinkScoreFiles

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
