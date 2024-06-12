import logging

from .lib import VariantFrame
from .lib import ScoringFileFrame, match_variants
from .lib import MatchResult, MatchResults
from .lib import PlinkScoreFiles


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

__version__ = "0.1.1"
