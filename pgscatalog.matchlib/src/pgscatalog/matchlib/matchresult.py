import collections.abc
import logging

import polars as pl

from pgscatalog.corelib import ZeroMatchesError

from ._plinkframe import PlinkFrames
from ._match.label import label_matches
from ._match.filter import filter_scores

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
    Container for MatchResult. Useful for making match logs and writing scoring files.

    >>> import tempfile, os, glob
    >>> from ._config import Config
    >>> from .variantframe import VariantFrame
    >>> from .scoringfileframe import ScoringFileFrame, match_variants
    >>> fout = tempfile.NamedTemporaryFile(delete=False)
    >>> target_path = Config.ROOT_DIR / "tests" / "good_match.pvar"
    >>> score_path =  Config.ROOT_DIR / "tests" / "good_match_scorefile.txt"
    >>> target = VariantFrame(target_path, dataset="goodmatch")
    >>> scorefile = ScoringFileFrame(score_path)
    >>> foutdir, splitfoutdir = tempfile.mkdtemp(), tempfile.mkdtemp()

    Using a context manager is import to prepare the data frames:

    >>> with target as target_df, scorefile as score_df:
    ...     results = match_variants(score_df=score_df, target_df=target_df, target=target)
    ...     _ = results.collect(outfile=fout.name)

    >>> with scorefile as score_df:
    ...     x = MatchResult.from_ipc(fout.name, dataset="goodmatch")
    ...     MatchResults(x).write_scorefiles(directory=foutdir, score_df=score_df)
    ...     MatchResults(x).write_scorefiles(directory=splitfoutdir, split=True, score_df=score_df)
    >>> MatchResults(x)  # doctest: +ELLIPSIS
    MatchResults([MatchResult(dataset=goodmatch, matchresult=None, ipc_path=...])

    By default, scoring files are written with multiple chromosomes per file:

    >>> combined_paths = glob.glob(foutdir + "/*ALL*")
    >>> combined_paths # doctest: +ELLIPSIS
    ['.../goodmatch_ALL_additive_0.scorefile.gz', '.../goodmatch_ALL_dominant_0.scorefile.gz', '.../goodmatch_ALL_recessive_0.scorefile.gz']
    >>> assert len(combined_paths) == 3

    Scoring files can be split. The input scoring file contains 20 unique
    chromosomes, with one additive + dominant effect file:

    >>> scorefiles = sorted(os.listdir(splitfoutdir))
    >>> scorefiles # doctest: +ELLIPSIS
    ['goodmatch_10_additive_0.scorefile.gz', 'goodmatch_11_additive_0.scorefile.gz', ...]
    >>> sum("dominant" in f for f in scorefiles)
    1
    >>> sum("recessive" in f for f in scorefiles)
    1
    >>> sum("additive" in f for f in scorefiles)
    19
    >>> assert len(scorefiles) == 21
    """

    def __init__(self, *elements):
        self._elements = list(elements)

        if len(datasets := set(x.dataset for x in elements)) > 1:
            raise ValueError(
                f"MatchResult elements must have same dataset, but found: {datasets!r}"
            )

        self.dataset = self._elements[0].dataset
        # a df composed of all match result elements
        self.df = pl.scan_ipc(x.ipc_path for x in self._elements)
        self._labelled = False

    def __len__(self):
        return len(self._elements)

    def __getitem__(self, item):
        return self._elements[item]

    def __repr__(self):
        return f"{type(self).__name__}({self._elements!r})"

    def label(self, **kwargs):
        """Label match candidates according to matching parameters

        kwargs control labelling parameters:

        **keep_first_match
        **remove_ambiguous
        **skip_flip
        **remove_multiallelic
        **filter_IDs
        """
        # if multiple match candidates are tied, keep the first (default: drop all)
        kwargs.setdefault("keep_first_match", False)
        # keep variants with ambiguous alleles, (e.g. A/T and G/C SNPs)
        kwargs.setdefault("remove_ambiguous", True)
        # consider matched variants that may be reported on the opposite strand
        kwargs.setdefault("skip_flip", False)
        # allow matching to multiallelic variants
        kwargs.setdefault("remove_multiallelic", True)
        # constrain variants to this list of IDs (ignored if empty list)
        kwargs.setdefault("filter_IDs", None)

        self._labelled = True
        return self.df.pipe(label_matches, kwargs)

    def write_scorefiles(self, directory, score_df, split=False, **kwargs):
        if not self._labelled:
            self.df = self.label(**kwargs)

        self.df, self.score_summary = filter_scores(
            scorefile=score_df,
            matches=self.df,
            min_overlap=kwargs.get("min_overlap", 0.75),
            dataset=self.dataset,
        )
        self.filtered = True

        if self.score_summary.is_empty():
            # can happen if min_overlap = 0
            raise ZeroMatchesError(
                "Error: no target variants match any variants in scoring files"
            )

        plink = PlinkFrames.from_matchresult(self.df)

        for frame in plink:
            frame.write(directory=directory, split=split, dataset=self.dataset)

    def full_variant_log(self):
        raise NotImplementedError
