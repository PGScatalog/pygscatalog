import logging
import importlib.metadata

from pgscatalog.match.lib import VariantFrame
from pgscatalog.match.lib import ScoringFileFrame, match_variants
from pgscatalog.match.lib import MatchResult, MatchResults
from pgscatalog.match.lib import PlinkScoreFiles


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

__version__ = importlib.metadata.version("pgscatalog.match")
