import os

import polars as pl

from pgscatalog.corelib import TargetVariants

from ._arrow import loose


class VariantFrame:
    """Similar to TargetVariants, but backed by the polars dataframe library

    Fast, supports more complicated things, but requires more resources (CPU/RAM)

    The context manager returns a polars LazyFrame:

    >>> path = Config.ROOT_DIR.parent / "pgscatalog.corelib" / "tests" / "hapnest.pvar"
    >>> x = VariantFrame(path)
    >>> with x as df:
    ...     df.collect().shape
    (100, 5)
    """

    def __init__(self, path, cleanup=True, tmpdir=None):
        self.variants = TargetVariants(path)
        self._cleanup = cleanup
        self._tmpdir = tmpdir
        self._loosed = False
        self.arrowpath = None

    def __enter__(self):
        # convert to arrow files
        self.arrowpath = loose(
            record_batches=self.variants.to_pa_recordbatch(),
            schema=self.variants.pa_schema(),
            tmpdir=self._tmpdir,
        )
        self._loosed = True
        # set up global string cache for categorical variables
        pl.enable_string_cache()
        return pl.scan_ipc(self.arrowpath.name)

    def __exit__(self, *args, **kwargs):
        if self._cleanup:
            os.unlink(self.arrowpath.name)
        self._loosed = False
        pl.disable_string_cache()
