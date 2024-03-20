import collections.abc
import logging

import polars as pl

from pgscatalog.core import ZeroMatchesError, MatchRateError

from ._plinkframe import PlinkFrames
from ._match.label import label_matches
from ._match.filter import filter_scores
from ._match.log import make_logs, check_log_count, make_summary_log

logger = logging.getLogger(__name__)


class MatchResult:
    """Represents variants in a scoring file matched against variants in a target genome

    When matching a scoring file, it's normal for matches to be composed of many
    :class:`MatchResult` objects. This is common if the target genome is split to
    have one chromosome per scoring file, and the container class
    :class:`MatchResults` provides some helpful methods for working with split data.

    >>> from ._config import Config
    >>> from .variantframe import VariantFrame
    >>> from .scoringfileframe import ScoringFileFrame, match_variants
    >>> target_path = Config.ROOT_DIR.parent / "pgscatalog.core" / "tests" / "data" / "hapnest.bim"
    >>> target = VariantFrame(target_path, dataset="hapnest")
    >>> score_path = Config.ROOT_DIR.parent / "pgscatalog.core" / "tests" / "data" / "combined.txt.gz"
    >>> scorefile = ScoringFileFrame(score_path)

    A :class:`MatchResult` can be instantiated with the lazyframe output of the match_variants function:

    >>> with target as target_df, scorefile as score_df:
    ...     match_variants(score_df=score_df, target_df=target_df, target=target)  # doctest: +ELLIPSIS
    MatchResult(dataset=hapnest, matchresult=[<LazyFrame...], ipc_path=None, df=None)

    A :class:`MatchResult` can also be saved to and loaded from Arrow IPC files:

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
        self.ipc_path = ipc_path
        # df: the dataframe backed by stable arrow files at ipc_path
        # after collect() is written to a file
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
        """Create an instance from an Arrow IPC file"""
        return cls(
            ipc_path=matchresults_ipc_path,
            dataset=dataset,
        )

    def __repr__(self):
        return f"{type(self).__name__}(dataset={self.dataset}, matchresult={self._matchresult!r}, ipc_path={self.ipc_path}, df={self.df!r})"


class MatchResults(collections.abc.Sequence):
    """
    Container for :class:`MatchResult`

    Useful for making matching logs and writing scoring files ready to be used by ``plink2 --score``

    >>> import tempfile, os, glob, pathlib
    >>> from ._config import Config
    >>> from .variantframe import VariantFrame
    >>> from .scoringfileframe import ScoringFileFrame, match_variants
    >>> fout = tempfile.NamedTemporaryFile(delete=False)
    >>> target_path = Config.ROOT_DIR / "tests" / "data" / "good_match.pvar"
    >>> score_path =  Config.ROOT_DIR / "tests" / "data" / "good_match_scorefile.txt"
    >>> target = VariantFrame(target_path, dataset="goodmatch")
    >>> scorefile = ScoringFileFrame(score_path)
    >>> foutdir, splitfoutdir = tempfile.mkdtemp(), tempfile.mkdtemp()

    Using a context manager is really important to prepare :class:`ScoringFileFrame` and :class:`VariantFrame` data frames:

    >>> with target as target_df, scorefile as score_df:
    ...     results = match_variants(score_df=score_df, target_df=target_df, target=target)
    ...     _ = results.collect(outfile=fout.name)

    These data frames are transparently backed by Arrow IPC files on disk.

    >>> with scorefile as score_df:
    ...     x = MatchResult.from_ipc(fout.name, dataset="goodmatch")
    ...     _ = MatchResults(x).write_scorefiles(directory=foutdir, score_df=score_df)
    ...     _ = MatchResults(x).write_scorefiles(directory=splitfoutdir, split=True, score_df=score_df)
    >>> MatchResults(x)  # doctest: +ELLIPSIS
    MatchResults([MatchResult(dataset=goodmatch, matchresult=None, ipc_path=...])

    By default, scoring files are written with multiple chromosomes per file:

    >>> combined_paths = sorted(glob.glob(foutdir + "/*ALL*"), key=lambda x: pathlib.Path(x).stem)
    >>> combined_paths # doctest: +ELLIPSIS
    ['.../goodmatch_ALL_additive_0.scorefile.gz', '.../goodmatch_ALL_dominant_0.scorefile.gz', '.../goodmatch_ALL_recessive_0.scorefile.gz']
    >>> assert len(combined_paths) == 3

    Scoring files can be split. The input scoring file contains 20 unique
    chromosomes, with one additive + dominant effect file (but one chromosome didn't match well):

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

    An important part of matching variants is reporting a log to see how well you're reproducing a PGS in the new target genomes:

    >>> with scorefile as score_df:
    ...     MatchResults(x).full_variant_log(score_df)  # +ELLIPSIS
    <LazyFrame [23 cols, {"row_nr": UInt64 … "dataset": Categorical(ordering='physical')}] at ...>
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
        if self.df.select("row_nr").collect().is_empty():
            raise ZeroMatchesError("No match candidates found for any scoring files")

        # a table containing up to one row per variant (the best possible match)
        self.filter_summary = None
        # a summary log containing match rates for variants
        self.summary_log = None
        # have match candidates in df been labelled?
        self._labelled = False
        # have match candidates in df been filtered?
        self._filtered = False
        # does the input scoring file count match the variant log count?
        self._log_OK = None

    def __len__(self):
        return len(self._elements)

    def __getitem__(self, item):
        return self._elements[item]

    def __repr__(self):
        return f"{type(self).__name__}({self._elements!r})"

    def label(
        self,
        keep_first_match=False,
        remove_ambiguous=True,
        skip_flip=False,
        remove_multiallelic=True,
        filter_IDs=None,
    ):
        """Label match candidates according to matching parameters

        kwargs control labelling parameters:

        * ``keep_first_match``: if best match candidates are tied, keep the first? (default: ```False``, drop all candidates for this variant)
        * ``remove_ambiguous``: Remove ambiguous alleles? (default: ``True``)
        * ``skip_flip``: Consider matched variants that may be reported on the opposite strand (default: ``False``)
        * ``remove_multiallelic`` remove multiallelic variants before matching (default: ``True``)
        * ``filter_IDs``: constrain variants to this list of IDs (default, don't constrain)
        """
        df = self.df.pipe(
            label_matches,
            keep_first_match=keep_first_match,
            remove_ambiguous=remove_ambiguous,
            skip_flip=skip_flip,
            remove_multiallelic=remove_multiallelic,
            filter_IDs=filter_IDs,
        )
        self._labelled = True
        self.df = df
        return self.df

    def filter(self, score_df, min_overlap=0.75, **kwargs):
        """Filter match candidates after labelling according to user parameters"""
        if not self._labelled:
            self.df = self.label(**kwargs)

        df, filter_summary = filter_scores(
            scorefile=score_df,
            matches=self.df,
            min_overlap=min_overlap,
            dataset=self.dataset,
        )
        self.filter_summary = filter_summary
        self.df = df
        self._filtered = True
        return self.df

    def write_scorefiles(
        self, directory, score_df, split=False, min_overlap=0.75, **kwargs
    ):
        """Write matches to a set of files ready for ``plink2 --score``

        Does some helpful stuff:

        * Labels match candidates
        * Filters match candidates based on labels and user configuration
        * Calculates match rates to see how well the PGS reproduces in the new target genomes
        * Generates a filtered variant log containing the best match candidate
        * Checks if the number of variants in the summary log matches the input scoring file
        * Sets up parallel score calculation (pivots data to wide column format)
        * Writes scores to a directory, splitting based on chromosome and effect type
        """
        if not self._labelled:
            _ = self.label(**kwargs)  # self.df gets updated

        if not self._filtered:
            _ = self.filter(score_df=score_df, min_overlap=min_overlap)

        if not all(x[1] for x in self.filter_summary.iter_rows()):
            logger.warning(f"{self.filter_summary}")
            raise MatchRateError(
                f"All scores fail to meet match threshold {min_overlap}"
            )

        # a summary log contains up to one variant (the best match) for each variant
        # in the scoring file
        self.summary_log = make_summary_log(
            match_candidates=self.df,
            dataset=self.dataset,
            filter_summary=self.filter_summary.lazy(),
            scorefile=score_df,
        )

        # double check log count vs scoring file variant count
        self._log_OK = check_log_count(scorefile=score_df, summary_log=self.summary_log)

        plink = PlinkFrames.from_matchresult(self.df)
        outfs = []
        for frame in plink:
            f = frame.write(directory=directory, split=split, dataset=self.dataset)
            outfs.append(f)

        # collect after joining in check_log_count (can't join df and lazy df)
        self.summary_log = self.summary_log.collect()
        return outfs

    def full_variant_log(self, score_df, **kwargs):
        """Generate a log for each variant in a scoring file

        Multiple match candidates may exist for each variant in the original file.
        Describe each variant (one variant per row) with match metadata
        """
        if not self._labelled:
            self.df = self.label(**kwargs)
            self._labelled = True

        if not self._filtered:
            self.df, self.filter_summary = filter_scores(
                scorefile=score_df,
                matches=self.df,
                min_overlap=kwargs.get("min_overlap", 0.75),
                dataset=self.dataset,
            )
            self._filtered = True

        return make_logs(
            scorefile=score_df, dataset=self.dataset, match_candidates=self.df
        )
