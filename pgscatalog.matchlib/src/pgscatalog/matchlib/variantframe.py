import os
import pathlib
import tempfile

from pgscatalog.corelib import TargetVariants

import pyarrow as pa


class VariantFrame:
    """Similar to TargetVariants, but backend runs with polars dataframe library

    Fast, supports more complicated things, but requires more resources (CPU/RAM)

    """

    def __init__(self, path):
        self._variants = TargetVariants(path)
        self._loosed = False
        self._arrowpath = None
        self._cleanup = None
        self._dir = None

    def loose(self, directory=None, cleanup=True):
        """Loose an arrow :) Stream compressed text files into temporary arrow files

        polars + arrow = very fast reading and processing

        >>> x = VariantFrame("/Users/bwingfield/Documents/projects/pygscatalog/pgscatalog.corelib/tests/hapnest.pvar")

        Cleaning up the temporary arrow files can be simplified by using a context manager:

        >>> import polars as pl
        >>> with x.loose() as arrow:
        ...     temppath = arrow._arrowpath.name
        ...     assert os.path.exists(temppath)
        ...     pl.read_ipc(temppath).shape
        (100, 5)
        >>> assert not os.path.exists(temppath)  # all cleaned up
        """
        self._cleanup = cleanup
        if directory is None:
            self._dir = tempfile.mkdtemp()
        else:
            self._dir = pathlib.Path(directory)

        self._arrowpath = tempfile.NamedTemporaryFile(dir=directory, delete=False)

        with pa.OSFile(self._arrowpath.name, "wb") as sink:
            with pa.RecordBatchFileWriter(
                sink=sink, schema=self._variants.pa_schema()
            ) as writer:
                for batch in self._variants.to_pa_recordbatch():
                    writer.write(batch)

        self._loosed = True
        return self

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._loosed = False
        if self._cleanup:
            os.unlink(self._arrowpath.name)
