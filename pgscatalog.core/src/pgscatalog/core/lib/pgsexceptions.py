""" This module defines a custom PGS exception hierarchy.

This module overrides sys.excepthook. Uncaught exceptions are intercepted,
and ``sys.exit()`` is called with a custom code. This can be helpful to debug
problems in other applications, because there are relatively few failure states.

The hierarchy:

- BasePGSException

  - MatchError

    - DuplicateMatchError

    - MatchRateError

    - ZeroMatchesError

    - MatchValueError

  - CombineError

    - BuildError

    - ScoreFormatError

    - LiftoverError

  - CatalogError

    - ScoreDownloadError

    - ScoreChecksumError

    - QueryError

    - InvalidAccessionError

  - SamplesheetError

    - GenomesNotFound

    - SamplesheetFormatError

"""
import sys
from types import MappingProxyType


class BasePGSException(Exception):
    """The base class from which all PGS errors must inherit.
    The purpose of this class is to simplify finding PGS exceptions and exiting python
    with a matching custom exit code."""

    def __init__(self, message):
        super().__init__(message)


class MatchError(BasePGSException):
    """The base class for errors that are raised during variant matching"""

    def __init__(self, message):
        super().__init__(message)


class DuplicateMatchError(MatchError):
    """Raised when a matched variant has been duplicated, so that a variant with the same ID
    would be split across two rows in an output scoring file.
    """


class MatchRateError(MatchError):
    """Raised when match rate is below match threshold for one or more scoring files"""


class ZeroMatchesError(MatchError):
    """Raised when zero matches are found for one or more scoring files.

    Distinct from MatchRateError because it's very common, and caused by bad input data or parameters."""


class MatchValueError(MatchError):
    """Raised when a match function receives inappropriate values.

    e.g., Multiple chromosomes detected in variant data but data is split per-chromosome"""


class CombineError(BasePGSException):
    """The base class for errors that are raised when combining scorefiles"""

    def __init__(self, message):
        super().__init__(message)


class BuildError(CombineError):
    """Raised when there's a problem with a scoring file genome build."""


class ScoreFormatError(CombineError):
    """Raised when there's a problem with a scoring file."""


class LiftoverError(CombineError):
    """Raised when liftover has failed to convert genomic coordinates well"""


class CatalogError(BasePGSException):
    """The base class for errors when querying or downloading from the PGS Catalog"""

    def __init__(self, message):
        super().__init__(message)


class ScoreDownloadError(CatalogError):
    """Raised when a scoring file can't be downloaded"""


class ScoreChecksumError(CatalogError):
    """Raised when a scoring file fails checksum validation"""


class QueryError(CatalogError):
    """Raised when the Catalog API doesn't return a valid response"""


class InvalidAccessionError(CatalogError):
    """Raised when an invalid term is used to query the Catalog"""


class SamplesheetError(BasePGSException):
    """The base class for errors related to samplesheet parsing"""

    def __init__(self, message):
        super().__init__(message)


class GenomesNotFound(SamplesheetError):
    """Raised when FileNotFound"""


class SamplesheetFormatError(SamplesheetError):
    """Raised when a samplesheet is badly formatted"""


class ExceptionExitCodeMap:
    """A read only map to get exit codes for custom exceptions"""

    # https://unix.stackexchange.com/a/604262
    _mapping = {
        ScoreDownloadError: 8,
        ScoreFormatError: 9,
        ScoreChecksumError: 10,
        QueryError: 11,
        InvalidAccessionError: 12,
        DuplicateMatchError: 13,
        MatchRateError: 14,
        ZeroMatchesError: 15,
        MatchValueError: 16,
        BuildError: 17,
        GenomesNotFound: 19,
        SamplesheetFormatError: 20,
    }

    code_map = MappingProxyType(_mapping)

    def __getitem__(self, exception_type):
        # if an exception can't be found in the map, return an error code (> 0) but default
        # max possible value 255
        return self.code_map.get(exception_type, 255)


def handle_uncaught_exception(exctype, value, trace):
    """Intercept BasePGSExceptions and trigger sys.exit with a custom code"""
    code_map = ExceptionExitCodeMap()
    oldHook(exctype, value, trace)
    if isinstance(value, BasePGSException):
        sys.exit(code_map[exctype])


sys.excepthook, oldHook = handle_uncaught_exception, sys.excepthook
