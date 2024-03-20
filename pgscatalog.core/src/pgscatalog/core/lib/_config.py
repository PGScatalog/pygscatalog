""" This module managers global package configuration"""

import importlib.metadata
import pathlib


class Config:
    """This class stores global package configuration as class attributes

    Most of the time you won't need to change anything.
    """

    # headers sent to the PGS Catalog API
    API_HEADER = {
        "user-agent": f"{__package__}/{importlib.metadata.version('pgscatalog.core')}"
    }
    # couldn't figure out a nicer way to get the root dir with a namespace package
    ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent.parent.parent.parent
    # when querying PGS Catalog API / downloading
    MAX_RETRIES = 5
    # disable HTTPS downloads from PGS Catalog
    FTP_EXCLUSIVE = False
    # the number of rows to read from a scoring file at a time
    BATCH_SIZE = 20000
    # the number of rows to read from a variant information file at a time
    TARGET_BATCH_SIZE = 500000
