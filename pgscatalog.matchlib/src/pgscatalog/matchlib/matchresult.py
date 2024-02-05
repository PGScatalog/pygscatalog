import collections.abc
import gzip
import io
import itertools
import logging
import pathlib

import polars as pl

from ._plinkframe import PlinkFrames
from ._match.label import label_matches

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

    >>> import tempfile, os
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
    >>> scorefiles = sorted(os.listdir(foutdir))
    >>> scorefiles # doctest: +ELLIPSIS
    ['goodmatch_10_additive_0.scorefile.gz', 'goodmatch_11_additive_0.scorefile.gz', ...]
    >>> sum("dominant" in f for f in scorefiles)
    1
    >>> sum("recessive" in f for f in scorefiles)
    1
    >>> sum("additive" in f for f in scorefiles)
    20
    >>> assert len(scorefiles) == 22
    >>> with scorefile as score_df:
    ...     MatchResults(x).summary_log(score_df=score_df)
    """

    def __init__(self, *elements):
        self._elements = list(elements)

        if len(datasets := set(x.dataset for x in elements)) > 1:
            raise ValueError(
                f"MatchResult elements must have same dataset, but found: {datasets!r}"
            )

        self.dataset = self._elements[0].dataset
        # a df composed of all match result elements
        self.df = pl.scan_ipc(x.ipc_path for x in self._elements).select(
            pl.col("*"), pl.col("match_type").cast(pl.Categorical)
        )
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
        # constrain variants to this list of IDs
        kwargs.setdefault("filter_IDs", [])

        self.df = self.df.pipe(label_matches, kwargs)
        self._labelled = True
        return self.df

    def write_scorefiles(self, directory, split=False, **kwargs):
        if not self._labelled:
            _ = self.label(**kwargs)

        # TODO: make plink frames from labelled and filtered data...
        plink = PlinkFrames.from_matchresult(self.df)
        # TODO: create summary log before writing - need to explode

        for frame in plink:
            frame.write(directory=directory, split=split, dataset=self.dataset)

    def summary_log(self, score_df, min_overlap=0.75):
        """ """
        # match_candidates = pl.concat(x.df for x in self._elements)
        # x, y = filter_scores(scorefile=score_df, dataset=self.dataset, matches=match_candidates, min_overlap=min_overlap)
        pass

    def variant_log(self):
        raise NotImplementedError


class PlinkScoreFiles(collections.abc.Sequence):
    """Represents a sequence of files written by MatchResults.write_scorefiles()

    Useful to combine split scoring files without re-computing matches
    """

    def __init__(self, *elements):
        self._elements = [pathlib.Path(x) for x in list(elements)]

    def __len__(self):
        return len(self._elements)

    def __getitem__(self, item):
        return self._elements[item]

    def __repr__(self):
        return f"{type(self).__name__}({self._elements!r})"

    def merge(self, directory):
        """Merge scoring files without recomputing matches

        Assumes a standard file naming system was used: dataset_chrom_effecttype_n

        >>> import tempfile, os, glob
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
        >>> foutdir = tempfile.mkdtemp()
        >>> MatchResults(x).write_scorefiles(directory=foutdir, split=True)  # doctest: +ELLIPSIS
        >>> plink_files = (pathlib.Path(foutdir) / x for x in os.listdir(foutdir))
        >>> psf = PlinkScoreFiles(*plink_files)
        >>> psf  # doctest: +ELLIPSIS
        PlinkScoreFiles([PosixPath('.../goodmatch_13_additive_0.scorefile.gz'), ...])
        >>> psf.merge(foutdir)
        >>> combined_paths = glob.glob(foutdir + "/*ALL*")
        >>> combined_paths # doctest: +ELLIPSIS
        ['.../goodmatch_ALL_additive_0.scorefile.gz', '.../goodmatch_ALL_dominant_0.scorefile.gz', '.../goodmatch_ALL_recessive_0.scorefile.gz']
        """
        dataset = self._elements[0].stem.split("_")[0]
        for x in self._elements:
            if dataset not in x.stem:
                raise ValueError(f"Invalid dataset: {dataset} and {x.stem}")

        def keyfunc(path):
            return path.stem.split("_")[2:]

        sorted_paths = sorted(self._elements, key=keyfunc)

        for k, g in itertools.groupby(sorted_paths, key=keyfunc):
            # multi-chrom -> ALL
            fout = "_".join([dataset, "ALL", *k]) + ".gz"
            paths = list(g)
            # infer_schema_length: read all columns as utf8 to simplify joining
            df = pl.concat(
                pl.read_csv(x, separator="\t", infer_schema_length=0) for x in paths
            ).fill_null(value="0")
            with gzip.open(pathlib.Path(directory) / fout, "wb") as gcsv:
                outf = io.TextIOWrapper(gcsv)
                df.write_csv(outf, separator="\t")
