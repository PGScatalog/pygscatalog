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
from typing import TYPE_CHECKING

import numpy as np
import zarr
import zarr.errors

from ..constants import (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    ZARR_COMPRESSOR,
    ZARR_VARIANT_CHUNK_SIZE,
)
from ..types import Pathish
from .zarrmodels import ZarrSampleMetadata, ZarrVariantMetadata

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    import numpy.typing as npt


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
        lengths = set(
            map(
                len,
                [self._chr_name, self._pos, self._refs, self._alts, gts],
            )
        )

        if len(lengths) != 1:
            raise ValueError(f"Input lists have different lengths: {lengths=}")

        # check that the genotypes look sensible
        if not np.all(
            np.isin(self._genotypes, [0, 1, MISSING_GENOTYPE_SENTINEL_VALUE])
        ):
            raise ValueError("Genotypes contain invalid values")

    def __len__(self) -> int:
        return len(self._chr_name)  # number of variants

    def __repr__(self):
        n_variants = len(self._pos)
        n_samples = len(self._samples)
        gt_shape = self._genotypes.shape
        return (
            f"<TargetVariants variants={n_variants} samples={n_samples} "
            f"genotypes_shape={gt_shape}"
        )

    @property
    def samples(self) -> list[str]:
        return self._samples

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
            **{
                "chr_name": self._chr_name,
                "chr_pos": self._pos,
                "ref": self._refs,
                "alts": self._alts,
                "sampleset": [self._sampleset] * len(self._chr_name),
                "filename": [str(self._target_path)] * len(self._chr_name),
                "variant_id": self.variant_ids,
            }
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
            gt_arr = zarr_group["genotypes"]
            gt_arr.append(self._genotypes)

        end_time = time.perf_counter()
        logger.info(
            f"{self._genotypes.shape} genotypes written to disk in "
            f"{end_time - start_time} seconds"
        )

    def _write_zarr_samples(self, zarr_group: zarr.Group) -> None:
        """Add sample metadata to the file group. A file must always have the same samples.

        Sample names in a sampleset may differ across files (e.g. a VCF might be split
        into batches of 100,000 samples).
        """
        sample_metadata: ZarrSampleMetadata = ZarrSampleMetadata(
            **{"samples": self.samples}
        )

        if "samples" not in zarr_group.attrs:
            logger.info(f"Adding {len(sample_metadata)} sample IDs to zarr attribute")
            zarr_group.attrs["samples"] = sample_metadata.model_dump()
        else:
            logger.info("Checking that sample IDs are consistent")
            existing_samples: ZarrSampleMetadata = ZarrSampleMetadata(
                **zarr_group.attrs["samples"]
            )
            if existing_samples.model_dump() != sample_metadata.model_dump():
                logger.critical(
                    f"Inconsistent sample IDs {self._target_path} (this should never happen)"
                )
                raise ValueError
            logger.info("Samples IDs are consistent")

    def _write_zarr_variants(self, zarr_group: zarr.Group) -> None:
        start_time = time.perf_counter()
        variant_metadata: ZarrVariantMetadata = self.variant_metadata

        if "variants" not in zarr_group.attrs:
            logger.info(
                f"Initialising variant metadata with {len(variant_metadata)} variants"
            )
            zarr_group.attrs["variants"] = variant_metadata.model_dump()
        else:
            logger.info("Updating existing variant metadata")
            existing_variants: ZarrVariantMetadata = ZarrVariantMetadata(
                **zarr_group.attrs["variants"]
            )
            logger.info(f"{len(existing_variants)} variants present in metadata")
            merged = existing_variants.merge(variant_metadata).model_dump()
            zarr_group.attrs["variants"] = merged
            logger.info(f"{len(merged)} variants present after metadata update")

        end_time = time.perf_counter()
        logger.info(
            f"{len(variant_metadata)} variants written to disk in "
            f"{end_time - start_time} seconds"
        )
