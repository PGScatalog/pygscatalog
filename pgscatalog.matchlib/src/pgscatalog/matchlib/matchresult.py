import collections.abc
import logging
import pathlib

import polars as pl

from .plinkframe import PlinkFrames
from ._match.plink import pivot_score

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

    def __repr__(self):
        return f"{type(self).__name__}(dataset={self.dataset}, matchresult={self._matchresult!r}, ipc_path={self.ipc_path}, df={self.df!r})"


class MatchResults(collections.abc.Sequence):
    """
    Container for MatchResult. Useful for making logs and writing scoring files.
    >>> import tempfile
    >>> from ._config import Config
    >>> from .variantframe import VariantFrame
    >>> from .scoringfileframe import ScoringFileFrame, match_variants
    >>> fout = tempfile.NamedTemporaryFile(delete=False)
    >>> target_path = Config.ROOT_DIR / "tests" / "good_match.pvar"
    >>> score_path =  Config.ROOT_DIR / "tests" / "good_match_scorefile.txt"
    >>> target = VariantFrame(target_path, dataset="goodmatch")
    >>> scorefile = ScoringFileFrame(score_path)

    >>> with target as target_df, scorefile as score_df:
    ...     results = match_variants(score_df=score_df, target_df=target_df, target=target)
    ...     _ = results.collect(outfile=fout.name)
    >>> x = MatchResult.from_ipc(fout.name, dataset="goodmatch")
    >>> MatchResults(x)  # doctest: +ELLIPSIS
    MatchResults([MatchResult(dataset=goodmatch, matchresult=None, ipc_path=...])
    >>> foutdir = tempfile.mkdtemp()
    >>> MatchResults(x).write_scorefiles(directory=foutdir)  # doctest: +ELLIPSIS
    >>> os.listdir(foutdir)
    ['goodmatch_ALL_additive_0.scorefile.gz', 'goodmatch_ALL_dominant_0.scorefile.gz', 'goodmatch_ALL_recessive_0.scorefile.gz']
    >>> assert len(os.listdir(foutdir)) == 3

    Scoring files can be split. The input scoring file contains 20 unique
    chromosomes, with one additive + dominant effect file:

    >>> foutdir = tempfile.mkdtemp()
    >>> MatchResults(x).write_scorefiles(directory=foutdir, split=True)  # doctest: +ELLIPSIS
    >>> x = sorted(os.listdir(foutdir))
    >>> x # doctest: +ELLIPSIS
    ['goodmatch_10_additive_0.scorefile.gz', 'goodmatch_11_additive_0.scorefile.gz', ...]
    >>> sum("dominant" in f for f in x)
    1
    >>> sum("recessive" in f for f in x)
    1
    >>> sum("additive" in f for f in x)
    20
    >>> assert len(x) == 22
    """

    def __init__(self, *elements):
        self._elements = list(elements)

        if len(datasets := set(x.dataset for x in elements)) > 1:
            raise ValueError(
                f"MatchResult elements must have same dataset, but found: {datasets!r}"
            )

        self.dataset = self._elements[0].dataset

    def __len__(self):
        return len(self._elements)

    def __getitem__(self, item):
        return self._elements[item]

    def __repr__(self):
        return f"{type(self).__name__}({self._elements!r})"

    def write_scorefiles(self, directory, split=False):
        plink = PlinkFrames.from_matchresult(self._elements)
        # TODO: create summary log before writing - need to explode

        for frame in plink:
            if split:
                dfs = frame.split()
                for chrom, df in dfs.items():
                    fout = (
                        pathlib.Path(directory)
                        / f"{self.dataset}_{chrom}_{str(frame.effect_type)}_{frame.n}.scorefile.gz"
                    )
                    df.pipe(pivot_score).write_csv(fout, separator="\t")
            else:
                chrom = "ALL"
                fout = (
                    pathlib.Path(directory)
                    / f"{self.dataset}_{chrom}_{str(frame.effect_type)}_{frame.n}.scorefile.gz"
                )
                frame.df.collect().pipe(pivot_score).write_csv(fout, separator="\t")

    def variant_log(self):
        raise NotImplementedError
