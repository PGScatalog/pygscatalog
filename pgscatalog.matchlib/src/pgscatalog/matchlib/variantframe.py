import os

from pgscatalog.corelib import TargetVariants

from ._arrow import loose


class VariantFrame:
    """Similar to TargetVariants, but backed by the polars dataframe library

    Fast, supports more complicated things, but requires more resources (CPU/RAM)
    """

    def __init__(self, path, cleanup=True, tmpdir=None):
        self._variants = TargetVariants(path)
        self._cleanup = cleanup
        self._tmpdir = tmpdir
        self._loosed = False
        self.arrowpath = None

    def __enter__(self):
        # convert to arrow files
        self.arrowpath = loose(
            record_batches=self._variants.to_pa_recordbatch(),
            schema=self._variants.pa_schema(),
            tmpdir=self._tmpdir,
        )
        self._loosed = True
        return self

    def __exit__(self, *args, **kwargs):
        if self._cleanup:
            os.unlink(self.arrowpath.name)
        self._loosed = False
