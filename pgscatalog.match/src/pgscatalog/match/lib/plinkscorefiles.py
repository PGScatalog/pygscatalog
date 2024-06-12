"""This module contains the PlinkScoreFiles class, which represents one or more
scoring files ready to be used with plink2 --score"""

import collections.abc
import gzip
import io
import itertools
import logging
import pathlib

import polars as pl

from pgscatalog.core import chrom_keyfunc, effect_type_keyfunc

logger = logging.getLogger(__name__)


class PlinkScoreFiles(collections.abc.Sequence):
    """Represents a sequence of scoring files written by :class:`MatchResults`"""

    def __init__(self, *elements):
        self._elements = [pathlib.Path(x) for x in sorted(list(elements))]

    def __len__(self):
        return len(self._elements)

    def __getitem__(self, item):
        return self._elements[item]

    def __repr__(self):
        return f"{type(self).__name__}({self._elements!r})"

    def merge(self, directory):
        """Merge scoring files without recomputing matches

        Assumes a standard file naming system was used: ``dataset_chrom_effecttype_n``

        >>> import tempfile, os, glob
        >>> from ._config import Config
        >>> from .variantframe import VariantFrame
        >>> from .scoringfileframe import ScoringFileFrame, match_variants
        >>> from .matchresult import MatchResult, MatchResults
        >>> fout = tempfile.NamedTemporaryFile(delete=False)
        >>> target_path = Config.ROOT_DIR / "tests" / "data" / "good_match.pvar"
        >>> score_path =  Config.ROOT_DIR / "tests" / "data" / "good_match_scorefile.txt"
        >>> target = VariantFrame(target_path, dataset="goodmatch")
        >>> scorefile = ScoringFileFrame(score_path)
        >>> with target as target_df, scorefile as score_df:
        ...     results = match_variants(score_df=score_df, target_df=target_df, target=target)
        ...     _ = results.collect(outfile=fout.name)
        >>> x = MatchResult.from_ipc(fout.name, dataset="goodmatch")
        >>> foutdir = tempfile.mkdtemp()
        >>> with scorefile as score_df:
        ...     _ = MatchResults(x).write_scorefiles(directory=foutdir, split=True, score_df=score_df)  # doctest: +ELLIPSIS
        >>> plink_files = (pathlib.Path(foutdir) / x for x in os.listdir(foutdir))
        >>> psf = PlinkScoreFiles(*plink_files)
        >>> psf  # doctest: +ELLIPSIS
        PlinkScoreFiles([PosixPath('.../goodmatch_1_additive_0.scorefile.gz'), ...])
        >>> psf.merge(foutdir)
        >>> combined_paths = sorted(glob.glob(foutdir + "/*ALL*"),  key=lambda x: pathlib.Path(x).stem)
        >>> len(combined_paths)
        3
        >>> combined_paths # doctest: +ELLIPSIS
        ['.../goodmatch_ALL_additive_0.scorefile.gz', '.../goodmatch_ALL_dominant_0.scorefile.gz', '.../goodmatch_ALL_recessive_0.scorefile.gz']

        """
        dataset = self._elements[0].stem.split("_")[0]
        for x in self._elements:
            if dataset not in x.stem:
                raise ValueError(f"Invalid dataset: {dataset} and {x.stem}")

        sorted_paths = sorted(self._elements, key=effect_type_keyfunc())

        for k, g in itertools.groupby(sorted_paths, key=effect_type_keyfunc()):
            logger.info(f"Writing combined scoring file {k}")

            # tidy up keys to create the output file name
            # keyfunc returns:
            # (('additive',), ('', 0.0, 'scorefile'))
            # need: [additive, 0]
            keys = (x for x in (itertools.chain(*k)) if x != "")
            keys = list(str(int(x)) if isinstance(x, float) else x for x in keys)
            keys.pop()  # drop scorefile

            # multi-chrom -> ALL
            fout = "_".join([dataset, "ALL", *keys]) + ".scorefile.gz"
            paths = sorted(list(g), key=chrom_keyfunc())

            # infer_schema_length: read all columns as utf8 to simplify joining
            dfs = (pl.read_csv(x, separator="\t", infer_schema_length=0) for x in paths)
            # diagonal concat is important to handle different column sets across dfs
            df = pl.concat(dfs, how="diagonal").fill_null(value="0")
            logger.info("Score files combined successfully")

            with gzip.open(
                pathlib.Path(directory) / fout, "wb", compresslevel=6
            ) as gcsv:
                logger.info(f"Writing out to {fout}")
                outf = io.TextIOWrapper(gcsv)
                df.write_csv(outf, separator="\t")
