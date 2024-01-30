import os

import polars as pl

from pgscatalog.corelib import NormalisedScoringFile

from ._arrow import loose


class ScoringFileFrame:
    """Like NormalisedScoringFile, but backed by the polars dataframe library

    Instantiated with a NormalisedScoringFile written to a file (i.e. the output of
    combine scorefiles application):

    >>> path = Config.ROOT_DIR.parent / "pgscatalog.corelib" / "tests" / "combined.txt.gz"
    >>> x = ScoringFileFrame(path)
    >>> x  # doctest: +ELLIPSIS
    ScoringFileFrame(NormalisedScoringFile('.../combined.txt.gz'))

    The context manager returns a polars dataframe with global string cache enabled:

    >>> with x as arrow:
    ...     assert os.path.exists(x.arrowpath.name)
    ...     arrow.collect().shape
    (154, 9)
    >>> assert not os.path.exists(x.arrowpath.name)
    """

    def __init__(self, path, cleanup=True, tmpdir=None):
        self.scoringfile = NormalisedScoringFile(path)
        self._cleanup = cleanup
        self._tmpdir = tmpdir
        self._loosed = False
        self.arrowpath = None

    def __enter__(self):
        # convert to arrow files
        self.arrowpath = loose(
            record_batches=self.scoringfile.to_pa_recordbatch(),
            schema=self.scoringfile.pa_schema(),
            tmpdir=self._tmpdir,
        )
        self._loosed = True
        pl.enable_string_cache()
        return pl.scan_ipc(self.arrowpath.name)

    def __exit__(self, *args, **kwargs):
        if self._cleanup:
            os.unlink(self.arrowpath.name)
        self._loosed = False
        pl.disable_string_cache()

    def __repr__(self):
        return f"{type(self).__name__}({repr(self.scoringfile)})"
