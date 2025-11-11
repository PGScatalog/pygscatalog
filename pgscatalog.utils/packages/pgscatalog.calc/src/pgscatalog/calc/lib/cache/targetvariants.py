"""
Provides the TargetVariants class for representing hard-called genotypes and variants

This module defines a lightweight container of variant information, including:

- chromosome name
- chromosome position
- reference allele
- alternate allele

Genotypes are stored in numpy unsigned 8-bit integer arrays (values between 0 - 255).
Valid genotype values include 0, 1, and a sentinel value to represent missing data.

The class exposes accessors for:

- genotypes: a 3D numpy array of shape (n_variants, n_samples, ploidy)
- samples: a list of sample identifiers
- variant_df: a polars dataframe which contains variant metadata

The module depends on:

- numpy for genotype matrix ops
- polars for building variant metadata tables quickly
"""

from __future__ import annotations

import logging
import time
from typing import TYPE_CHECKING, cast

import numpy as np
import zarr
import zarr.errors

from pgscatalog.calc.lib.constants import (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    ZARR_COMPRESSOR,
    ZARR_VARIANT_CHUNK_SIZE,
)

from .zarrmodels import ZarrSampleMetadata, ZarrVariantMetadata

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from collections.abc import Sized

    import numpy.typing as npt

    from pgscatalog.calc.lib.types import Pathish


class TargetVariants:
    def __init__(
        self,
        chr_name: list[str],
        pos: list[int],
        refs: list[str | None],
        alts: list[list[str] | None],
        gts: list[npt.NDArray[np.uint8]],
        samples: list[str],
        target_path: Pathish,
        sampleset: str,
    ):
        self._chr_name = chr_name
        self._pos = pos
        self._refs = refs
        self._alts = alts
        self._genotypes = np.stack(gts)
        self._samples = samples
        self._target_path = target_path
        self._sampleset = sampleset

        # check that all input lists have the same length
        to_check: list[Sized] = [self._chr_name, self._pos, self._refs, self._alts, gts]
        lengths = {len(x) for x in to_check}

        if len(lengths) != 1:
            raise ValueError(f"Input lists have different lengths: {lengths=}")

        # check that the genotypes look sensible
        if not np.all(
            np.isin(self._genotypes, [0, 1, MISSING_GENOTYPE_SENTINEL_VALUE])
        ):
            raise ValueError("Genotypes contain invalid values")

    def __len__(self) -> int:
        return len(self._chr_name)  # number of variants

    def __repr__(self) -> str:
        n_variants = len(self._pos)
        n_samples = len(self._samples)
        gt_shape = self._genotypes.shape
        return (
            f"<TargetVariants variants={n_variants} samples={n_samples} "
            f"genotypes_shape={gt_shape}"
        )

    @property
    def samples(self) -> list[str]:
        return list(self._samples)

    @property
    def genotypes(self) -> npt.NDArray[np.uint8]:
        return self._genotypes

    @property
    def variant_ids(self) -> list[str]:
        def _alt_to_str(alt: list[str] | None) -> str:
            if alt is None:
                return ""

            return "/".join(alt)

        ids = []
        for chrom, pos, ref, alt in zip(
            self._chr_name, self._pos, self._refs, self._alts, strict=True
        ):
            ids.append(f"{chrom}:{pos}:{ref or ''}:{_alt_to_str(alt)}")

        return ids

    @property
    def variant_metadata(self) -> ZarrVariantMetadata:
        """Convert variant metadata to a dict"""
        logger.info("Converting TargetVariants to polars dataframe")

        return ZarrVariantMetadata(
            chr_name=list(self._chr_name),
            chr_pos=list(self._pos),
            ref=list(self._refs),
            alts=list(self._alts),
            variant_id=list(self.variant_ids),
        )

    def write_zarr(self, zarr_group: zarr.Group) -> None:
        """Write TargetVariants to a zarr group

        Sample IDs, variant metadata, and a genotype array is written to the zarr group

        The group must be at a file level in the hierarchy
        """
        self._write_zarr_samples(zarr_group)
        self._write_zarr_variants(zarr_group)
        self._write_zarr_array(zarr_group)

    def _write_zarr_array(self, zarr_group: zarr.Group) -> None:
        """Store a 3D genotype array in the genotypes group"""
        start_time = time.perf_counter()
        try:
            logger.info("Creating zarr array")
            zarr_chunks = (
                ZARR_VARIANT_CHUNK_SIZE,
                self.genotypes.shape[1],
                self.genotypes.shape[2],
            )

            zarr_group.create_array(
                chunks=zarr_chunks,
                data=self._genotypes,
                fill_value=MISSING_GENOTYPE_SENTINEL_VALUE,
                name="genotypes",
                compressors=ZARR_COMPRESSOR,
            )
            logger.info(f"Zarr array {self.genotypes.shape=} created")
        except zarr.errors.ContainsArrayError:
            # array has already been initialised
            logger.info("Array already exists, appending")
            gt_arr = cast("zarr.Array", zarr_group["genotypes"])
            gt_arr.append(self._genotypes)

        end_time = time.perf_counter()
        logger.info(
            f"{self._genotypes.shape} genotypes written to disk in "
            f"{end_time - start_time} seconds"
        )

    def _write_zarr_samples(self, zarr_group: zarr.Group) -> None:
        """Add sample metadata to the file group. A file must always have the same
        samples.

        Sample names in a sampleset may differ across files (e.g. a VCF might be split
        into batches of 100,000 samples).
        """
        sample_metadata: ZarrSampleMetadata = ZarrSampleMetadata.model_validate(
            self.samples
        )

        if "samples" not in zarr_group.attrs:
            logger.info(f"Adding {len(sample_metadata)} sample IDs to zarr attribute")
            zarr_group.attrs["samples"] = sample_metadata.model_dump()
        else:
            logger.info("Checking that sample IDs are consistent")
            existing_samples: ZarrSampleMetadata = ZarrSampleMetadata.model_validate(
                zarr_group.attrs["samples"]
            )

            if existing_samples.model_dump() != sample_metadata.model_dump():
                logger.critical(
                    f"Inconsistent sample IDs {self._target_path} "
                    f"(this should never happen)"
                )
                raise ValueError
            logger.info("Samples IDs are consistent")

    def _write_zarr_variants(self, zarr_group: zarr.Group) -> None:
        """Store variant metadata as a set of 1D arrays in the meta group

        This means that variant metadata supports compression and the numpy arrays
        can be read into polars

        See https://github.com/zarr-developers/zarr-specs/discussions/365
        """
        meta_arrays = self.variant_metadata.to_numpy()
        zarr_chunks = (ZARR_VARIANT_CHUNK_SIZE,)
        meta_root = zarr.open_group(
            store=zarr_group.store, path=f"{zarr_group.path}/meta"
        )
        logger.info(f"Writing variant metadata to zarr group {meta_root.info}")

        # mypy can't work out what to do when iterating over a TypedDict
        for metadata, array in meta_arrays.items():
            if self.genotypes.shape[0] != array.shape[0]:
                raise ValueError(
                    "1D variant metadata array must have same number of rows as "
                    "3D genotype array"
                )

            start_time = time.perf_counter()
            try:
                logger.info(f"Creating variant {metadata=} zarr array")
                meta_root.create_array(
                    chunks=zarr_chunks,
                    data=array,
                    name=metadata,
                    compressors=ZARR_COMPRESSOR,
                )
                logger.info(f"Zarr array {array.shape=} created")
            except zarr.errors.ContainsArrayError:
                # array has already been initialised
                logger.info("Array already exists, appending")
                metadata_arr = cast("zarr.Array", meta_root[metadata])
                metadata_arr.append(array)

            end_time = time.perf_counter()
            logger.info(
                f"{metadata} {array.shape=} written to disk in "
                f"{end_time - start_time} seconds"
            )


def add_missing_positions_to_lists(
    *,
    chroms: list[str],
    positions: list[int],
    ref_alleles: list[str | None],
    alt_alleles: list[list[str] | None],
    hard_calls: list[npt.NDArray[np.uint8]],
    scoring_file_regions: list[tuple[str, int]],
    seen_positions: set[tuple[str, int]],
    n_samples: int,
) -> None:
    """Mutates lists in place!"""
    # represent missing genotypes with a sentinel value of the same shape
    missing_gts = np.full(
        shape=(n_samples, 2),
        fill_value=MISSING_GENOTYPE_SENTINEL_VALUE,
        dtype=np.uint8,
    )

    missing_positions = set(scoring_file_regions) - seen_positions
    for chr_name, chr_pos in missing_positions:
        # important to prevent future cache misses
        chroms.append(chr_name)
        positions.append(chr_pos)
        ref_alleles.append(None)
        alt_alleles.append(None)
        hard_calls.append(missing_gts)
