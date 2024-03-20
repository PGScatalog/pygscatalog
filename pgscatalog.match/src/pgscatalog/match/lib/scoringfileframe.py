import logging
import os
import shutil

import polars as pl

from pgscatalog.core import NormalisedScoringFile

from ._arrow import loose
from ._match.preprocess import complement_valid_alleles
from ._match.match import get_all_matches
from .matchresult import MatchResult

logger = logging.getLogger(__name__)


def match_variants(score_df, target_df, target):
    """Get all match candidates for a VariantFrame dataframe and ScoringFileFrame dataframe

    Returns a :class:`MatchResult`
    """
    matches = get_all_matches(scorefile=score_df, target=target_df)
    return MatchResult(dataset=target.dataset, matchresult=matches)


class ScoringFileFrame:
    """Like :class:`pgscatalog.core.NormalisedScoringFile`, but backed by the polars dataframe library

    Instantiated with a :class:`pgscatalog.core.NormalisedScoringFile` written to
    a file. This is a long format/melted CSV file containing normalised variant data
    (i.e. the output of combine scorefiles application):

    >>> from ._config import Config
    >>> path = Config.ROOT_DIR.parent / "pgscatalog.core" / "tests" / "data" / "combined.txt.gz"
    >>> x = ScoringFileFrame(path)
    >>> x  # doctest: +ELLIPSIS
    ScoringFileFrame(NormalisedScoringFile('.../combined.txt.gz'))

    Using a context manager is important to prepare a polars dataframe:

    >>> with x as arrow:
    ...     assert os.path.exists(x.arrowpath.name)
    ...     arrow.collect().shape
    (154, 11)
    >>> assert not os.path.exists(x.arrowpath.name)  # all cleaned up

    >>> from .variantframe import VariantFrame
    >>> path = Config.ROOT_DIR.parent / "pgscatalog.core" / "tests" / "data" / "hapnest.bim"
    >>> target = VariantFrame(path, dataset="hapnest")
    >>> with target as target_df, x as score_df:
    ...     match_variants(score_df=score_df, target_df=target_df, target=target)  # doctest: +ELLIPSIS
    MatchResult(dataset=hapnest, matchresult=[<LazyFrame [16 cols...
    """

    def __init__(self, path, chrom=None, cleanup=True, tmpdir=None):
        self.scoringfile = NormalisedScoringFile(path)
        # used for filtering the scoring file if the target contains a single chromosome
        self.chrom = chrom
        self._cleanup = cleanup
        self._tmpdir = tmpdir
        self._loosed = False
        self.arrowpath = None

    def __enter__(self):
        """Use a context manager to create arrow IPC files and return a lazyframe"""
        # prevent nested context managers creating multiple arrow files
        pl.enable_string_cache()

        if not self._loosed:
            logger.debug(f"Converting {self!r} to feather format")
            self.arrowpath = loose(
                record_batches=self.scoringfile.to_pa_recordbatch(),
                schema=self.scoringfile.pa_schema(),
                tmpdir=self._tmpdir,
            )
            self._loosed = True
            logger.debug(f"{self!r} feather conversion complete")

        if self.arrowpath is not None:
            score_df = (
                pl.scan_ipc(self.arrowpath.name)
                .pipe(
                    complement_valid_alleles,
                    flip_cols=["effect_allele", "other_allele"],
                )
                .with_columns(
                    [
                        pl.col("chr_name").cast(pl.Categorical),
                        pl.col("accession").cast(pl.Categorical),
                        pl.col("effect_allele").cast(pl.Categorical),
                        pl.col("other_allele").cast(pl.Categorical),
                        pl.col("effect_allele_FLIP").cast(pl.Categorical),
                        pl.col("other_allele_FLIP").cast(pl.Categorical),
                    ]
                )
            )
        else:
            raise ValueError("No arrow IPC file")

        if self.chrom is not None:
            # add filter to query plan
            logger.debug(f"Filtering scoring file to chromosome {self.chrom}")
            score_df = score_df.filter(pl.col("chr_name") == self.chrom)

        return score_df

    def __exit__(self, *args, **kwargs):
        """Optionally clean up the arrow IPC files"""
        if self._cleanup:
            os.unlink(self.arrowpath.name)
            self._loosed = False
        pl.disable_string_cache()

    def __repr__(self):
        return f"{type(self).__name__}({repr(self.scoringfile)})"

    def save_ipc(self, destination):
        """Save the dataframe prepared by the context manager to an Arrow IPC file

        Useful because the context manager will clean up the IPC files while exiting.

        This method allows data to be persisted."""
        if not self._loosed:
            raise ValueError(
                "Can't save IPC because it doesn't exist."
                "Try calling inside a with: statement"
            )

        shutil.copyfile(self.arrowpath, destination)
