from __future__ import annotations

import logging
from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np
import polars as pl

from ._bgen_probabilities import (
    phased_probabilities_to_hard_calls,
    unphased_probabilities_to_hard_calls,
)
from .constants import (
    BGEN_PHASED_N_COLS,
    BGEN_UNPHASED_N_COLS,
    MISSING_GENOTYPE_SENTINEL_VALUE,
    ZARR_PLOIDY,
)
from .genomefiletypes import GenomeFileType

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from collections.abc import Iterable

    import dask.array as da
    import numpy.typing as npt


@dataclass(slots=True)
class TargetVariant:
    """A small set of fields from a single row in a VCF"""

    chr_name: str
    chr_pos: int
    gts: npt.NDArray[np.uint8] | None = None
    probs: npt.NDArray[np.float64] | None = None
    is_phased: bool | None = None
    ref: str | None = None
    alts: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        if self.gts is not None and self.gts.shape[-1] > 2:
            raise ValueError("Array shape is not haploid or diploid")

        # remove any padding if it exists
        self.chr_name = self.chr_name.lstrip("0")

    @property
    def variant_id(self) -> str:
        ref = self.ref if self.ref else ""
        alts = "/".join(self.alts) if self.alts else ""
        return f"{self.chr_name}:{self.chr_pos}:{ref}:{alts}"


class TargetVariants:
    """A collection of TargetVariant

    Provides methods to convert to a columnar format (arrow / parquet)
    """

    def __init__(
        self,
        variants: Iterable[TargetVariant],
        samples: list[str],
        filename: str,
        sampleset: str,
        filetype: GenomeFileType,
    ):
        self._variants = list(variants)
        self._samples = list(set(samples))  # thinking about duplicated ids
        self.filetype = filetype

        try:
            gts_n_samples: int = list(
                {v.gts.shape[0] for v in self._variants if v.gts is not None}
            )[0]
        except IndexError:
            # skip this check when all variants are missing genotypes
            pass
        else:
            n_samples: int = len(self._samples)
            if gts_n_samples != n_samples:
                raise ValueError(
                    f"Number of samples in genotypes ({gts_n_samples=} doesn't match "
                    f"provided sample length {n_samples=})"
                )

        self._filename = filename  # the file that the variants were collected from
        self._sampleset = sampleset  # a human friendly label

    def __len__(self) -> int:
        return len(self.variants)

    @property
    def variants(self) -> list[TargetVariant]:
        return self._variants

    @property
    def samples(self) -> list[str]:
        return self._samples

    @property
    def variant_ids(self) -> pl.Series:
        return pl.Series([str(x.variant_id) for x in self.variants], dtype=pl.String)

    @cached_property
    def genotypes(self) -> npt.NDArray[np.uint8]:
        """
        Return genotypes as a stacked 3D matrix.

        Shape: (variants, samples, ploidy)
        - Axis 0: Variants
        - Axis 1: Samples
        - Axis 2: Ploidy

        This mirrors the layout of `to_polars()`.
        """
        # dummy value for None list elements
        null_gt_calls = np.full(
            (len(self.samples), ZARR_PLOIDY),
            fill_value=MISSING_GENOTYPE_SENTINEL_VALUE,
            dtype=np.uint8,
        )
        cleaned_arrays = [
            (x.gts if x.gts is not None else null_gt_calls) for x in self._variants
        ]
        return np.stack(cleaned_arrays, axis=0)

    def probs_to_hard_calls(self) -> da.Array:
        is_phased = set(x.is_phased for x in self._variants)
        if True in is_phased and False in is_phased:
            raise ValueError("Variants having mixed phase")

        # None phashing values are missing gts, must be filled by dispatched functions
        is_any_phased = True in is_phased
        n_samples = len(self.samples)

        logger.info(f"{self.filetype} {is_phased=}")
        match (self.filetype, is_any_phased):
            case (GenomeFileType.BGEN, True):
                # phasing: { None, True } or { True }
                missing = np.full(
                    (n_samples, BGEN_PHASED_N_COLS),
                    fill_value=MISSING_GENOTYPE_SENTINEL_VALUE,
                    dtype=np.uint8,
                )
                probs = [
                    x.probs if x.probs is not None else missing for x in self._variants
                ]
                hard_calls = phased_probabilities_to_hard_calls(
                    probabilities=probs, n_samples=n_samples
                )
            case (GenomeFileType.BGEN, False):
                # phasing: { None, False}, { False }, or { None } (all variants
                # missing gt)
                missing = np.full(
                    (n_samples, BGEN_UNPHASED_N_COLS),
                    fill_value=MISSING_GENOTYPE_SENTINEL_VALUE,
                    dtype=np.uint8,
                )
                probs = [
                    x.probs if x.probs is not None else missing for x in self._variants
                ]
                hard_calls = unphased_probabilities_to_hard_calls(
                    probabilities=probs, n_samples=n_samples
                )
            case (GenomeFileType.VCF, _):
                raise NotImplementedError
            case _:
                raise NotImplementedError

        return hard_calls

    @cached_property
    def variant_df(self) -> pl.DataFrame:
        """Return a polars dataframe from an iterator of targetvariants

        columnar format, so let's convert our iterable into a dictionary
        where keys = column names and values = lists (rows)
        """
        logger.info("Converting TargetVariants to polars dataframe")
        fields = ["chr_name", "chr_pos", "ref", "alts"]
        merged: dict[str, list] = {field: [] for field in fields}

        for variant in self.variants:
            for field in fields:
                merged[field].append(getattr(variant, field))

        df = pl.DataFrame(merged).with_columns(
            pl.col("chr_name").cast(pl.String),
            pl.col("chr_pos").cast(pl.UInt32),
            pl.col("ref").cast(pl.String),
            pl.col("alts").cast(pl.List(pl.String)),
            sampleset=pl.lit(self._sampleset),
            filename=pl.lit(self._filename),
            variant_id=self.variant_ids,
        )

        logger.info("Polars dataframe ready")
        return df
