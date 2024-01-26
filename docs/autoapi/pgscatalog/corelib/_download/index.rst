:orphan:

:py:mod:`pgscatalog.corelib._download`
======================================

.. py:module:: pgscatalog.corelib._download

.. autoapi-nested-parse::

   This module contains functions that simplify downloading data from the PGS Catalog.
   HTTPS download is preferred, with FTP fallback available.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pgscatalog.corelib._download.ftp_fallback
   pgscatalog.corelib._download.https_download
   pgscatalog.corelib._download.score_download_failed



Attributes
~~~~~~~~~~

.. autoapisummary::

   pgscatalog.corelib._download.logger


.. py:function:: ftp_fallback(retry_state)

   Try downloading from PGS Catalog using FTP like it's 1991 instead of HTTPS.
   This function should be used as a callback function from a tenacity.retry decorator.

   Downloading over FTP is less reliable but a helpful fallback option.
   It should never be called directly, only as a callback function from tenacity.retry:

   >>> with tempfile.TemporaryDirectory() as d:
   ...     out_path = "PGS000001.txt.gz"
   ...     kwargs = {"url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000001/ScoringFiles/PGS000001.txt.gz", "directory": d, "out_path": out_path}
   ...     retry_state = tenacity.RetryCallState(retry_object=None, fn=None, args=None, kwargs=kwargs)
   ...     ftp_fallback(retry_state)
   ...     assert (pathlib.Path(d) / pathlib.Path(out_path)).exists()



.. py:function:: https_download(*, url, out_path, directory, overwrite)

   Download a file from the PGS Catalog over HTTPS, with automatic retries and
   waiting. md5 checksums are automatically validated.


.. py:function:: score_download_failed(retry_state)

   This function raises a ``ScoreChecksumError`` or ``ScoreDownloadError``
   because every attempt has failed. This function should be used as a callback
   function from a tenacity.retry decorator


.. py:data:: logger

   

