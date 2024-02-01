import logging

import polars as pl

logger = logging.getLogger(__name__)


class MatchResult:
    """This class represents variant match results

    >>> from ._config import Config
    >>> from .variantframe import VariantFrame
    >>> from .scoringfileframe import ScoringFileFrame, match_variants
    >>> target_path = Config.ROOT_DIR.parent / "pgscatalog.corelib" / "tests" / "hapnest.bim"
    >>> target = VariantFrame(target_path, dataset="hapnest")
    >>> score_path = Config.ROOT_DIR.parent / "pgscatalog.corelib" / "tests" / "combined.txt.gz"
    >>> scorefile = ScoringFileFrame(score_path)

    MatchResults can be instantiated with the lazyframe output of match_variants:

    >>> with target as target_df, scorefile as score_df:
    ...     match_variants(score_df=score_df, target_df=target_df, target=target)  # doctest: +ELLIPSIS
    MatchResult(dataset=hapnest, matchresult=[<LazyFrame...], ipc_path=None, df=None)

    MatchResults can also be loaded from IPC files:

    >>> import tempfile
    >>> fout = tempfile.NamedTemporaryFile(delete=False)
    >>> with target as target_df, scorefile as score_df:
    ...     results = match_variants(score_df=score_df, target_df=target_df, target=target)
    ...     _ = results.collect(outfile=fout.name)
    >>> x = MatchResult.from_ipc(fout.name, dataset="hapnest")
    >>> x  # doctest: +ELLIPSIS
    MatchResult(dataset=hapnest, matchresult=None, ipc_path=..., df=<LazyFrame...>)
    >>> x + x # doctest: +ELLIPSIS
    MatchResult(dataset=hapnest, matchresult=[<LazyFrame...], ipc_path=None, df=None)
    """

    def __init__(self, dataset, matchresult=None, ipc_path=None):
        # a dataset ID
        self.dataset = dataset
        # df: the dataframe backed by stable arrow files at ipc_path
        # after collect() is written to a file
        self.ipc_path = ipc_path
        self.df = None

        # a lazyframe of match results, backed by arrow files that may get cleaned up
        self._matchresult = matchresult

        match (self.ipc_path, matchresult):
            case (_, None):
                # init from results already written to a file
                self.df = pl.scan_ipc(self.ipc_path)
            case (None, _):
                # init from a lazy frame result from get_all_matches()
                # needs to be collected() later and optionally written to an IPC file
                pass
            case (None, None) | (_, _):
                raise ValueError

    def collect(self, outfile=None):
        """Compute match results and optionally save to file"""
        if self._matchresult is None:
            raise ValueError("Can't collect, missing matchresult")

        if outfile is not None:
            pl.concat(pl.collect_all(self._matchresult)).write_ipc(outfile)
            self.df = pl.scan_ipc(outfile)
        else:
            # call .lazy() to prevent eager queries later
            self.df = pl.concat(pl.collect_all(self._matchresult)).lazy()

        return self.df

    @classmethod
    def from_ipc(cls, matchresults_ipc_path, dataset):
        return cls(
            ipc_path=matchresults_ipc_path,
            dataset=dataset,
        )

    def __add__(self, other):
        if not isinstance(other, MatchResult):
            return NotImplemented

        if self.dataset != other.dataset:
            logger.warning(f"Dataset mismatch: {self.dataset}, {other.dataset}")
            return NotImplemented

        if self.df is None or other.df is None:
            logger.warning("Stable dataframe is missing. Did you .collect()?")
            return NotImplemented

        return MatchResult(
            matchresult=[pl.concat([self.df, other.df]).lazy()], dataset=self.dataset
        )

    def __repr__(self):
        return f"{type(self).__name__}(dataset={self.dataset}, matchresult={self._matchresult!r}, ipc_path={self.ipc_path}, df={self.df!r})"
