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
from typing import TYPE_CHECKING

import numpy as np
import polars as pl

from ..constants import MISSING_GENOTYPE_SENTINEL_VALUE
from ..types import Pathish

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    import numpy.typing as npt


class TargetVariants:
    def __init__(
        self,
        chr_name: list[str],
        pos: list[int],
        refs: list[str],
        alts: list[str],
        gts: list[npt.NDArray[np.uint8]],
        samples: list[str],
        target_path: Pathish,
        sampleset: str
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

    def metadata_to_df(self, sampleset: str) -> pl.DataFrame:
        """ Convert variant metadata to a polars dataframe """
        logger.info("Converting TargetVariants to polars dataframe")

        merged = {}
        merged["chr_name"] = self._chr_name
        merged["chr_pos"] = self._pos
        merged["ref"] = self._refs
        merged["alts"] = self._alts

        df = pl.DataFrame(merged).with_columns(
            pl.col("chr_name").cast(pl.String),
            pl.col("chr_pos").cast(pl.UInt32),
            pl.col("ref").cast(pl.String),
            pl.col("alts").cast(pl.List(pl.String)),
            sampleset=pl.lit(self._sampleset),
            filename=pl.lit(str(self._target_path)),
            variant_id=pl.concat_str(
                [pl.col("chr_name"), pl.col("chr_pos"), pl.col("ref"), pl.col("alts")],
                separator=":",
                ignore_nulls=True,
            ),
        )

        logger.info("Polars dataframe ready")
        return df
