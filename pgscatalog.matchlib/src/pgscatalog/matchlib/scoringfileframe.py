import os

from pgscatalog.corelib import NormalisedScoringFile

from ._arrow import loose


class ScoringFileFrame:
    """Like NormalisedScoringFile, but backed by the polars dataframe library

    Instantiated with a NormalisedScoringFile written to a file (i.e. the output of
    combine scorefiles application):

    >>> import polars as pl
    >>> path = Config.ROOT_DIR.parent / "pgscatalog.corelib" / "tests" / "combined.txt.gz"
    >>> x = ScoringFileFrame(path)
    >>> x  # doctest: +ELLIPSIS
    ScoringFileFrame(NormalisedScoringFile('.../combined.txt.gz'))
    >>> with x as arrow:
    ...     assert os.path.exists(arrow.arrowpath.name)
    ...     pl.read_ipc(arrow.arrowpath.name).shape
    (154, 9)
    >>> assert not os.path.exists(arrow.arrowpath.name)
    """

    def __init__(self, path, cleanup=True, tmpdir=None):
        self._scoringfile = NormalisedScoringFile(path)
        self._cleanup = cleanup
        self._tmpdir = tmpdir
        self._loosed = False
        self.arrowpath = None

    def __enter__(self):
        # convert to arrow files
        self.arrowpath = loose(
            record_batches=self._scoringfile.to_pa_recordbatch(),
            schema=self._scoringfile.pa_schema(),
            tmpdir=self._tmpdir,
        )
        self._loosed = True
        return self

    def __exit__(self, *args, **kwargs):
        if self._cleanup:
            os.unlink(self.arrowpath.name)
        self._loosed = False

    def __repr__(self):
        return f"{type(self).__name__}({repr(self._scoringfile)})"
