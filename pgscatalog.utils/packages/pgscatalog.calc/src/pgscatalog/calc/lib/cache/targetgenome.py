from __future__ import annotations

import logging
import os
import pathlib
from functools import cached_property
from typing import TYPE_CHECKING

import duckdb
import polars as pl
import zarr
import zarr.errors
import zarr.storage

from ._genomefilehandlers import GenomeFileHandler, get_file_handler
from .genomefiletypes import GenomeFileType
from .targetvariants import TargetVariants
from .zarrmodels import get_position_df

if TYPE_CHECKING:
    from collections.abc import Sequence

    from ..types import Pathish


logger = logging.getLogger(__name__)


class TargetGenome:
    def __init__(
        self,
        *,
        target_path: Pathish,
        cache_dir: Pathish,
        sampleset: str,
        chrom: str | None = None,
        threads: int = 2,
        sample_file: Pathish | None = None,
    ) -> None:
        self._target_path = target_path
        self._cache_dir = pathlib.Path(cache_dir).resolve()
        self._sampleset = sampleset
        self._chrom = chrom
        self._threads = threads
        self._handler: GenomeFileHandler = get_file_handler(
            path=target_path,
            sample_file=sample_file,
            cache_dir=self._cache_dir,
            sampleset=self._sampleset,
        )

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"target_path={os.fspath(self._target_path)!r}, "
            f"chrom={self._chrom!r}, "
            f"cache_dir={os.fspath(self._cache_dir) if self._cache_dir else None!r})"
        )

    @property
    def sampleset(self) -> str:
        """A human label for a (set of) target genomes, e.g. UKBiobank"""
        return self._sampleset

    @property
    def filename(self) -> str:
        return pathlib.Path(self._target_path).name

    @property
    def filetype(self) -> GenomeFileType:
        return self._handler.file_type

    @property
    def cache_dir(self) -> pathlib.Path:
        return pathlib.Path(self._cache_dir)

    @cached_property
    def samples(self) -> list[str]:
        return self._handler.samples

    @property
    def target_path(self) -> pathlib.Path:
        return pathlib.Path(self._target_path)

    @property
    def _zarr_archive_name(self):
        """The name of the directory containing the zarr array root"""
        return self.cache_dir / "genotypes.zarr"

    @property
    def _zarr_group_path(self):
        """The group hierarchy:

        / (root of zarr array
        └── <sampleset> zarr group
            └── <filename> zarr group
                └── genotypes (n_variants, n_samples, PLOIDY) uint8 array
        """
        return f"{self.sampleset}/{pathlib.Path(self.filename).name}"

    @property
    def _zarr_store(self) -> zarr.storage.StoreLike:
        """A directory on the local filesystem contains the zarr array"""
        return zarr.storage.LocalStore(root=self._zarr_archive_name, read_only=False)

    @property
    def zarr_group(self) -> zarr.Group:
        return zarr.group(
            store=self._zarr_store,
            path=self._zarr_group_path,
            overwrite=False,
            zarr_version=3,
            zarr_format=3,
        )

    @property
    def chrom(self) -> str | None:
        return self._chrom

    @property
    def cached_positions(self) -> pl.DataFrame:
        """A dataframe of positions previously queried

        Contains two columns: chr_name and chr_pos
        """
        return get_position_df(self.zarr_group)

    def cache_variants(self, positions: Sequence[tuple[str, int]]) -> None:
        """Query an indexed target genome and store variants in a zarr array"""
        logger.info(f"{len(positions)} positions requested")

        # no support for X / XY / Y / patches yet
        # drop them here but keep them in scoring files to correctly calculate overlap
        with duckdb.connect() as conn:
            cached = self.cached_positions
            logger.info(f"{cached.shape[0]} variants cached")

            query = pl.DataFrame(
                positions, schema={"chr_name": pl.String, "chr_pos": pl.Int64}
            )
            logger.info(f"{query.shape[0]} rows to query")

            missing: list[tuple[str, int]] = conn.sql("""
                WITH cached AS (SELECT * FROM cached),
                    queries AS (SELECT * FROM query),
                    missing AS (SELECT * FROM queries
                        ANTI JOIN cached 
                        USING (chr_name, chr_pos))
                SELECT chr_name, chr_pos
                FROM missing
                WHERE chr_name IN (
                    SELECT CAST(i AS VARCHAR)
                    FROM UNNEST(range(1, 23)) AS t(i) 
                )
                ORDER BY (chr_name, chr_pos);
            """).fetchall()
            logger.info(f"{len(missing)} rows missing from cache")

        if len(missing) == 0:
            logger.info("No positions to query")
            return

        variants: TargetVariants = self._handler.query_variants(positions=missing)

        # variants may be empty with very small batches or sparse files
        if variants:
            variants.write_zarr(zarr_group=self.zarr_group)
        else:
            logger.warning("No variants to write. This is weird!")

        logger.info(f"{self} finished caching")
