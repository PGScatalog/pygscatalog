""" This module contains global corelib configuration stored as class attributes."""

import importlib.metadata
import pathlib


class Config:
    API_HEADER = {
        "user-agent": f"{__package__}/{importlib.metadata.version(__package__)}"
    }
    # couldn't figure out a nicer way to get the root dir with a namespace package
    ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent.parent.parent
    # when querying PGS Catalog API / downloading
    MAX_RETRIES = 5
    # disable HTTPS downloads from PGS Catalog
    FTP_EXCLUSIVE = False
    # the number of rows to read from a scoring file at a time
    BATCH_SIZE = 20000
