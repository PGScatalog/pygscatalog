"""
This configuration class stores settings as class attributes

Some zarr settings need to be configurable at runtime based on CLI arguments.

It's more reliable to import this config instance and update clas attributes compared
with importing the constants module and overriding values there.
"""

from __future__ import annotations

import zarr.codecs


class Config:
    ZARR_VARIANT_CHUNK_SIZE: int = 5000
    ZARR_COMPRESSOR = zarr.codecs.ZstdCodec(level=5)


config = Config()
