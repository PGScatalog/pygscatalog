import logging
import os
import shutil

import polars as pl

from pgscatalog.core import TargetVariants

from ._arrow import loose
from ._match.preprocess import filter_target, annotate_multiallelic

logger = logging.getLogger(__name__)


class VariantFrame:
    """Similar to :class:`pgscatalog.core.TargetVariants`, but backed by the polars dataframe library

    Fast, supports more complicated things, but requires more resources.

    The context manager returns a polars LazyFrame:

    >>> from ._config import Config
    >>> path = Config.ROOT_DIR.parent / "pgscatalog.core" / "tests" / "data" / "hapnest.bim"
    >>> x = VariantFrame(path, dataset="hapnest")
    >>> with x as df:
    ...     df.collect().shape
    (101, 6)
    >>> x  # doctest: +ELLIPSIS
    VariantFrame(path='.../hapnest.bim', dataset='hapnest', chrom=None, cleanup=True, tmpdir=None)

    The :class:`VariantFrame` contains a :class:`pgscatalog.core.TargetVariants` object:

    >>> x.variants  # doctest: +ELLIPSIS
    TargetVariants(path='.../hapnest.bim')
    """

    def __init__(self, path, dataset, chrom=None, cleanup=True, tmpdir=None):
        self.variants = TargetVariants(path)
        # variant frames need an ID column (dataset)
        self.dataset = dataset
        self.chrom = chrom
        self._cleanup = cleanup
        self._tmpdir = tmpdir
        self._loosed = False
        self.arrowpath = None

    def __enter__(self):
        if not self._loosed:
            # convert to arrow files
            self.arrowpath = loose(
                record_batches=self.variants.to_pa_recordbatch(),
                schema=self.variants.pa_schema(),
                tmpdir=self._tmpdir,
            )
            self._loosed = True

        # set up global string cache for categorical variables
        pl.enable_string_cache()

        target_df = (
            pl.scan_ipc(self.arrowpath.name)
            .pipe(filter_target)
            .pipe(annotate_multiallelic)
            .with_columns(
                [
                    pl.col("#CHROM").cast(pl.Categorical),
                    pl.col("REF").cast(pl.Categorical),
                    pl.col("ALT").cast(pl.Categorical),
                ]
            )
        )

        return target_df

    def __exit__(self, *args, **kwargs):
        if self._cleanup:
            os.unlink(self.arrowpath.name)
            self._loosed = False
        pl.disable_string_cache()

    def __repr__(self):
        return (
            f"{type(self).__name__}(path={repr(self.variants.path)}, dataset={repr(self.dataset)}, "
            f"chrom={repr(self.chrom)}, cleanup={repr(self._cleanup)}, "
            f"tmpdir={repr(self._tmpdir)})"
        )

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
