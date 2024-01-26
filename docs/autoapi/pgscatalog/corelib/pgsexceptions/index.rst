:orphan:

:py:mod:`pgscatalog.corelib.pgsexceptions`
==========================================

.. py:module:: pgscatalog.corelib.pgsexceptions

.. autoapi-nested-parse::

   This module defines a custom PGS exception hierarchy.

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



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   pgscatalog.corelib.pgsexceptions.ExceptionExitCodeMap



Functions
~~~~~~~~~

.. autoapisummary::

   pgscatalog.corelib.pgsexceptions.handle_uncaught_exception



.. py:exception:: BasePGSException(message)


   Bases: :py:obj:`Exception`

   The base class from which all PGS errors must inherit.
   The purpose of this class is to simplify finding PGS exceptions and exiting python
   with a matching custom exit code.


.. py:exception:: BuildError(message)


   Bases: :py:obj:`CombineError`

   Raised when there's a problem with a scoring file genome build.


.. py:exception:: CatalogError(message)


   Bases: :py:obj:`BasePGSException`

   The base class for errors when querying or downloading from the PGS Catalog


.. py:exception:: CombineError(message)


   Bases: :py:obj:`BasePGSException`

   The base class for errors that are raised when combining scorefiles


.. py:exception:: DuplicateMatchError(message)


   Bases: :py:obj:`MatchError`

   Raised when a matched variant has been duplicated, so that a variant with the same ID
   would be split across two rows in an output scoring file.


.. py:exception:: GenomesNotFound(message)


   Bases: :py:obj:`SamplesheetError`

   Raised when FileNotFound


.. py:exception:: InvalidAccessionError(message)


   Bases: :py:obj:`CatalogError`

   Raised when an invalid term is used to query the Catalog


.. py:exception:: LiftoverError(message)


   Bases: :py:obj:`CombineError`

   Raised when liftover has failed to convert genomic coordinates well


.. py:exception:: MatchError(message)


   Bases: :py:obj:`BasePGSException`

   The base class for errors that are raised during variant matching


.. py:exception:: MatchRateError(message)


   Bases: :py:obj:`MatchError`

   Raised when match rate is below match threshold for one or more scoring files


.. py:exception:: MatchValueError(message)


   Bases: :py:obj:`MatchError`

   Raised when a match function receives inappropriate values.

   e.g., Multiple chromosomes detected in variant data but data is split per-chromosome


.. py:exception:: QueryError(message)


   Bases: :py:obj:`CatalogError`

   Raised when the Catalog API doesn't return a valid response


.. py:exception:: SamplesheetError(message)


   Bases: :py:obj:`BasePGSException`

   The base class for errors related to samplesheet parsing


.. py:exception:: SamplesheetFormatError(message)


   Bases: :py:obj:`SamplesheetError`

   Raised when a samplesheet is badly formatted


.. py:exception:: ScoreChecksumError(message)


   Bases: :py:obj:`CatalogError`

   Raised when a scoring file fails checksum validation


.. py:exception:: ScoreDownloadError(message)


   Bases: :py:obj:`CatalogError`

   Raised when a scoring file can't be downloaded


.. py:exception:: ScoreFormatError(message)


   Bases: :py:obj:`CombineError`

   Raised when there's a problem with a scoring file.


.. py:exception:: ZeroMatchesError(message)


   Bases: :py:obj:`MatchError`

   Raised when zero matches are found for one or more scoring files.

   Distinct from MatchRateError because it's very common, and caused by bad input data or parameters.


.. py:class:: ExceptionExitCodeMap


   A read only map to get exit codes for custom exceptions

   .. py:attribute:: code_map

      


.. py:function:: handle_uncaught_exception(exctype, value, trace)

   Intercept BasePGSExceptions and trigger sys.exit with a custom code 


