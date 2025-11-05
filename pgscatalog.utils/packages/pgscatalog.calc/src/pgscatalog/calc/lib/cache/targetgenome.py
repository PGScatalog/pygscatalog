from __future__ import annotations

import logging
import os
import pathlib
import time
from functools import cached_property
from typing import TYPE_CHECKING

import duckdb
import numpy as np
import tenacity
import zarr
import zarr.errors
from filelock import BaseFileLock, FileLock

from ._genomefilehandlers import GenomeFileHandler, get_file_handler
from .genomefiletypes import GenomeFileType
from ..constants import (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    ZARR_COMPRESSOR,
    ZARR_MAX_N_VARIANTS,
    ZARR_PLOIDY,
    ZARR_VARIANT_CHUNK_SIZE,
)
from .targetvariants import TargetVariants

if TYPE_CHECKING:
    from collections.abc import Sequence

    import numpy.typing as npt
    import polars as pl

    from .types import Pathish


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
        """A directory that points to hive partitioned parquet files"""
        return pathlib.Path(self._cache_dir)

    @cached_property
    def samples(self) -> list[str]:
        return self._handler.samples

    @property
    def target_path(self) -> pathlib.Path:
        return pathlib.Path(self._target_path)

    @property
    def zarr_group(self) -> zarr.Group:
        filename = pathlib.Path(self.filename).name
        group_path = f"{self.sampleset}/{filename}"
        return zarr.group(store=self.cache_dir / "gts", path=group_path)

    @property
    def zarr_group_lock(self) -> BaseFileLock:
        lock_path = self.cache_dir / "gts" / self.zarr_group.path / ".lock"
        return FileLock(lock_path)

    @property
    def variants_db_path(self) -> pathlib.Path:
        return self.cache_dir / "variants.db"

    @property
    def chrom(self) -> str | None:
        return self._chrom

    @property
    @tenacity.retry(
        retry=tenacity.retry_if_exception_type(duckdb.IOException),
        wait=tenacity.wait_random_exponential(multiplier=1, max=60),
    )
    def cached_positions(self) -> frozenset[tuple[str, int]]:
        """A list of positions previously queried from a VCF"""
        cached_positions = []
        logger.info("Attempting to connect to database to check cache")

        try:
            with duckdb.connect(self.variants_db_path) as conn:
                logger.info("Connected to database")
                # filtering by filename is important to handle split files properly
                query = conn.sql(f"""
                    SELECT DISTINCT chr_name, chr_pos
                    FROM variants_table
                    WHERE filename = '{self.filename}'
                """)

                if self.chrom is not None:
                    query.filter(f"chr_name = {self.chrom}")

                cached_positions = query.fetchall()
        except duckdb.IOException as e:
            if "Cannot open file" in str(e):
                cached_positions = []
            else:
                raise duckdb.IOException("Busy") from e
        except duckdb.CatalogException:
            # db doesn't exist but that's OK
            cached_positions = []

        logger.info("Released database connection")
        return frozenset(cached_positions)

    def cache_variants(self, positions: Sequence[tuple[str, int]]) -> None:
        logger.info(f"{len(positions)} positions requested")

        # no support for X / XY / Y / patches yet
        # drop them here but keep them in scoring files to correctly calculate overlap
        valid_chromosomes: set[str] = {str(x) for x in range(1, 23)}
        positions_to_query: frozenset[tuple[str, int]] = frozenset(
            filter(lambda p: p[0] in valid_chromosomes, positions)
        )
        logger.info(
            f"Dropped {len(positions) - len(positions_to_query)} positions from "
            f"unsupported chromosomes (X/XY/MT/patch)"
        )

        logger.info(f"Caching {len(positions_to_query)} unique positions {self.chrom=}")
        sorted_positions: list[tuple[str, int]] = sorted(
            positions_to_query - self.cached_positions,
            key=lambda x: (int(x[0]), x[1]),
        )
        logger.info(
            f"Fetching {len(sorted_positions)} positions after checking local cache"
        )

        if len(sorted_positions) == 0:
            logger.info("No positions to query")
            return

        variants: TargetVariants = self._handler.query_variants(
            positions=sorted_positions
        )

        # variants may be empty with very small batches or sparse files
        if variants:
            self._write_variants(variants=variants)
        else:
            logger.warning("No variants to write. This is weird!")

        logger.info(f"{self} finished caching")

    def _write_variants(self, variants: TargetVariants) -> None:
        gt_indices = write_variants_to_db(
            database_path=self.variants_db_path,
            target_variants=variants,
            sampleset=self._sampleset,
        )
        write_genotypes_to_zarr(
            zarr_group=self.zarr_group,
            target_variants=variants,
            gt_indices=gt_indices,
            lock=self.zarr_group_lock,
        )
        write_metadata_to_zarr(
            cache_dir=self.cache_dir,
            sampleset=self.sampleset,
            samples=self.samples,
            filename=self.filename,
        )


def write_metadata_to_zarr(
    cache_dir: Pathish, sampleset: str, samples: list[str], filename: Pathish
) -> None:
    logger.info(f"Adding n_sample metadata to group {sampleset}")
    sampleset_group = zarr.group(
        store=pathlib.Path(cache_dir) / "gts", path=f"{sampleset}"
    )
    sample_list = sampleset_group.attrs.get("samples", None)

    if sample_list is None:
        sampleset_group.attrs["samples"] = samples
    else:
        if not isinstance(sample_list, list) or not all(
            isinstance(s, str) for s in sample_list
        ):
            raise TypeError("samples attribute must be a list of strings")
        bad_samples = set(sample_list).symmetric_difference(set(samples))
        if bad_samples:
            logger.critical(f"Bad samples detected: {bad_samples=}")
            raise ValueError(
                f"{len(sample_list)=} in metadata mismatches {len(samples)=}"
                f" in {filename=} "
            )


def init_genotype_array(zarr_group: zarr.Group, n_samples: int) -> zarr.Array:
    try:
        # rows, cols, z-axis
        zarr_shape = (
            ZARR_MAX_N_VARIANTS,
            n_samples,
            ZARR_PLOIDY,
        )
        # chunks must not overlap across processes, so set chunk size to batch size
        # (batch size is equivalent to length of target variants)
        zarr_chunks = (
            ZARR_VARIANT_CHUNK_SIZE,
            n_samples,
            ZARR_PLOIDY,
        )

        gt_arr: zarr.Array = zarr_group.create_array(
            shape=zarr_shape,
            dtype=np.uint8,
            chunks=zarr_chunks,
            fill_value=MISSING_GENOTYPE_SENTINEL_VALUE,
            name="genotypes",
            compressors=ZARR_COMPRESSOR,
        )
    except zarr.errors.ContainsArrayError:
        if not isinstance(zarr_group["genotypes"], zarr.Array):
            raise TypeError("Expected a Zarr Array >:(") from None
        # array has already been initialised
        gt_arr = zarr_group["genotypes"]

    return gt_arr


def write_genotypes_to_zarr(
    zarr_group: zarr.Group,
    target_variants: TargetVariants,
    gt_indices: npt.NDArray[np.int64],
    lock: BaseFileLock,
) -> None:
    """Vectorised function that writes genotypes to a zarr store

    A lock is used here because it's too difficult to align work over chunks
    """
    logger.info("Acquiring lock to write genotypes")
    with lock:
        logger.info("Genotype write lock acquired")
        start_time = time.perf_counter()
        zarr_array = init_genotype_array(
            zarr_group=zarr_group, n_samples=len(target_variants.samples)
        )
        zarr_array[gt_indices] = target_variants.genotypes
        end_time = time.perf_counter()

    logger.info(
        f"{target_variants.genotypes.shape} genotypes written to disk in "
        f"{end_time - start_time} seconds, lock released"
    )


@tenacity.retry(
    retry=tenacity.retry_if_exception_type(duckdb.IOException),
    wait=tenacity.wait_random_exponential(multiplier=1, max=60),
)
def write_variants_to_db(
    database_path: Pathish, target_variants: TargetVariants, sampleset: str
) -> npt.NDArray[np.int64]:
    start_time = time.perf_counter()
    sequence_name = f"sequence_{sampleset}"
    logger.info("Attempting to connect to database")
    with duckdb.connect(str(database_path)) as conn:
        logger.info("Connected to database")
        # each sampleset has its own sequence generator to get a monotonically
        # increasing integer. this is useful for setting and getting array indices
        # (genotypes)
        conn.sql(f"CREATE SEQUENCE IF NOT EXISTS {sequence_name};")
        conn.sql("""
        CREATE TABLE IF NOT EXISTS variants_table (
            variant_id TEXT NOT NULL,
            chr_name TEXT NOT NULL,
            chr_pos UINTEGER NOT NULL,
            ref TEXT,
            alts TEXT[],
            sampleset TEXT NOT NULL,
            filename TEXT NOT NULL,
            geno_index UINTEGER NOT NULL,
            PRIMARY KEY (variant_id, filename)
        );
        """)
        _df: pl.DataFrame = target_variants.metadata_to_df(sampleset=sampleset)
        try:
            conn.begin()
            indices = conn.execute(f"""
                    INSERT INTO variants_table
                    SELECT
                        variant_id,
                        chr_name,
                        chr_pos,
                        ref,
                        alts,
                        sampleset,
                        parse_filename(filename) AS filename,
                        nextval(\'{sequence_name}\') - 1 AS geno_index
                    FROM _df
                    RETURNING geno_index;
                    """).fetchall()
            conn.commit()
        except Exception:
            conn.rollback()
            raise

    end_time = time.perf_counter()
    logger.info(
        f"Inserted {_df.shape[0]} variants to database in "
        f"{(end_time - start_time):.2f} seconds"
    )
    logger.info("db connection released")
    return np.array(indices, dtype=np.int64).flatten()
