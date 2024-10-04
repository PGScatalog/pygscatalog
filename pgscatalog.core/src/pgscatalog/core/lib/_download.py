"""This module contains functions that simplify downloading data from the PGS Catalog.
HTTPS download is preferred, with FTP fallback available.
"""

import hashlib
import logging
import os
import pathlib
import tempfile
import urllib

import tenacity
import httpx

from pgscatalog.core.lib.pgsexceptions import ScoreDownloadError, ScoreChecksumError
from pgscatalog.core.lib._config import Config

logger = logging.getLogger(__name__)


def score_download_failed(retry_state):
    """This function raises a ``ScoreChecksumError`` or ``ScoreDownloadError``
    because every attempt has failed. This function should be used as a callback
    function from a tenacity.retry decorator"""
    try:
        retry_state.outcome.result()
    except ScoreChecksumError as e:
        raise ScoreChecksumError("All checksum retries failed") from e
    except Exception as download_exc:
        raise ScoreDownloadError("All download retries failed") from download_exc
    finally:
        logger.critical(f"Score download failed after all retries: {retry_state!r}")


@tenacity.retry(
    stop=tenacity.stop_after_attempt(Config.MAX_RETRIES),
    retry=tenacity.retry_if_exception_type((ScoreDownloadError, ScoreChecksumError)),
    retry_error_callback=score_download_failed,
    wait=tenacity.wait_fixed(3) + tenacity.wait_random(0, 2),
)
def ftp_fallback(retry_state):
    """Try downloading from PGS Catalog using FTP like it's 1991 instead of HTTPS.
    This function should be used as a callback function from a tenacity.retry decorator.

    Downloading over FTP is less reliable but a helpful fallback option.
    It should never be called directly, only as a callback function from tenacity.retry:

    >>> with tempfile.TemporaryDirectory() as d:
    ...     out_path = "PGS000001.txt.gz"
    ...     kwargs = {"url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000001/ScoringFiles/PGS000001.txt.gz", "directory": d, "out_path": out_path}
    ...     retry_state = tenacity.RetryCallState(retry_object=None, fn=None, args=None, kwargs=kwargs)
    ...     ftp_fallback(retry_state)
    ...     assert (pathlib.Path(d) / pathlib.Path(out_path)).exists()

    """
    url = retry_state.kwargs.get("url")
    directory = retry_state.kwargs.get("directory")

    ftp_url = url.replace("https://", "ftp://")
    checksum_url = ftp_url + ".md5"

    fn = pathlib.Path(retry_state.kwargs.get("out_path")).name
    out_path = pathlib.Path(directory) / pathlib.Path(fn)
    md5 = hashlib.md5()

    try:
        with (
            tempfile.NamedTemporaryFile(dir=directory, delete=False) as score_f,
            tempfile.NamedTemporaryFile(dir=directory, delete=True) as checksum_f,
        ):
            urllib.request.urlretrieve(ftp_url, score_f.name)
            urllib.request.urlretrieve(checksum_url, checksum_f.name)

            md5.update(score_f.read())

            if (checksum := md5.hexdigest()) != (
                remote := checksum_f.read().decode().split()[0]
            ):
                raise ScoreChecksumError(
                    f"Local checksum {checksum} doesn't match remote {remote}"
                )
    except urllib.error.URLError as download_exc:
        raise ScoreDownloadError(f"Can't download {ftp_url} over FTP") from download_exc
    else:
        # no exceptions thrown, move the temporary file to the final output path
        os.rename(score_f.name, out_path)
        logger.info(f"FTP download OK, {out_path} checksum validation passed")


@tenacity.retry(
    stop=tenacity.stop_after_attempt(Config.MAX_RETRIES),
    retry=tenacity.retry_if_exception_type((ScoreDownloadError, ScoreChecksumError)),
    retry_error_callback=ftp_fallback,
    wait=tenacity.wait_fixed(3) + tenacity.wait_random(0, 2),
)
def https_download(*, url, out_path, directory, overwrite):
    """Download a file from the PGS Catalog over HTTPS, with automatic retries and
    waiting. md5 checksums are automatically validated."""
    try:
        if Config.FTP_EXCLUSIVE:
            logger.warning("HTTPS downloads disabled by Config.FTP_EXCLUSIVE")
            https_download.retry.wait = tenacity.wait_none()
            https_download.retry.stop = tenacity.stop.stop_after_attempt(1)
            raise ScoreDownloadError("HTTPS disabled")

        if out_path.exists() and not overwrite:
            raise FileExistsError(f"{out_path} already exists")

        checksum_path = url + ".md5"
        checksum = httpx.get(checksum_path, headers=Config.API_HEADER).text
        md5 = hashlib.md5()

        with tempfile.NamedTemporaryFile(dir=directory, delete=False) as f:
            with httpx.stream("GET", url, headers=Config.API_HEADER) as r:
                for data in r.iter_bytes():
                    f.write(data)
                    md5.update(data)

            if (calc := md5.hexdigest()) != (remote := checksum.split()[0]):
                # will attempt to download again (see decorator)
                raise ScoreChecksumError(
                    f"Calculated checksum {calc} doesn't match {remote}"
                )
    except httpx.UnsupportedProtocol as protocol_exc:
        raise ValueError(f"Can't download a local file: {url!r}") from protocol_exc
    except httpx.RequestError as download_exc:
        raise ScoreDownloadError("HTTPS download failed") from download_exc
    else:
        # no exceptions thrown, move the temporary file to the final output path
        os.rename(f.name, out_path)
        logger.info(f"HTTPS download OK, {out_path} checksum validation passed")
