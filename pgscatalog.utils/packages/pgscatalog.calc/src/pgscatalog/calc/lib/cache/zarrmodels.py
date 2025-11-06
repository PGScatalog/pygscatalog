"""
Pydantic models for loading and saving data from zarr attributes (metadata)

These models help to read and write structured data about variants and samples
This data is helpful when working with the genotype array (variants = row names,
samples = column names).
"""

from __future__ import annotations

import logging
from copy import deepcopy
from typing import TYPE_CHECKING, Annotated, TypedDict

import polars as pl
import zarr
import zarr.errors
from pydantic import AfterValidator, BaseModel, PositiveInt, RootModel, model_validator

from calc.lib.constants import NUMPY_STRING_DTYPE

if TYPE_CHECKING:
    import numpy as np
    import numpy.typing as npt

logger = logging.getLogger(__name__)


class ZarrVariantMetadata(BaseModel):
    """A dataframe-y model, suitable for ingesting with pandas / polars / databases

    Useful for saving / loading data about variants into a group-level zarr attribute.

    Each target genome file will have its own variant metadata.
    """

    chr_name: list[str]
    chr_pos: list[PositiveInt]
    ref: Annotated[list[str | None], AfterValidator(is_valid_allele)]
    alts: Annotated[list[list[str] | None], AfterValidator(is_valid_allele)]
    variant_id: list[str]

    def __len__(self):
        return len(self.chr_name)

    @model_validator(mode="after")
    def check_even_length(self) -> ZarrVariantMetadata:
        lengths = [len(v) for v in self.__dict__.values()]
        if len(set(lengths)) > 1:
            raise ValueError(
                f"All list fields must have the same length, got {lengths}"
            )
        return self

    def merge(self, other: ZarrVariantMetadata) -> ZarrVariantMetadata:
        d1, d2 = self.model_dump(), other.model_dump()
        # deepcopy: any future changes won't affect original objects
        merged = {k: deepcopy(d1[k]) + deepcopy(d2[k]) for k in d1}
        return ZarrVariantMetadata(**merged)

    def to_df(self) -> pl.DataFrame:
        return pl.DataFrame(self.model_dump())

    def to_numpy(self) -> VariantArrays:
        """Convert the dataframe into 1D arrays

        Handles converting strings to a consistent fixed width dtype
        """
        df = self.to_df()

        # store strings as fixed-width: 255 unicode characters max (same as plink)
        chr_name = df["chr_name"].to_numpy().astype(NUMPY_STRING_DTYPE)
        chr_pos = df["chr_pos"].to_numpy()
        # ref may be null, numpy will unhelpfully convert None -> "None"
        ref = df["ref"].fill_null("").to_numpy().astype(NUMPY_STRING_DTYPE)
        # squash list column with / characters
        alt = (
            df["alts"]
            .list.join("/")
            .fill_null("")
            .to_numpy()
            .astype(NUMPY_STRING_DTYPE)
        )
        variant_id = df["variant_id"].to_numpy().astype(NUMPY_STRING_DTYPE)

        return {
            "chr_name": chr_name,
            "chr_pos": chr_pos,
            "ref": ref,
            "alts": alt,
            "variant_id": variant_id,
        }


class VariantArrays(TypedDict):
    chr_name: npt.NDArray[np.str_]
    chr_pos: npt.NDArray[np.int64]
    ref: npt.NDArray[np.str_]
    alts: npt.NDArray[np.str_]
    variant_id: npt.NDArray[np.str_]


def is_valid_allele(alleles: list[list[str] | None] | list[str | None]):
    valid_alleles = {"A", "C", "T", "G"}
    for allele in alleles:
        if allele is not None and not set(allele).issubset(valid_alleles):
            raise ValueError(f"Invalid allele: {allele}")
    return alleles


def get_position_df(zarr_group: zarr.Group) -> pl.DataFrame:
    """Get variants from the "meta" zarr group array"""
    try:
        meta_root = zarr.open_group(
            path=f"{zarr_group.path}/meta", store=zarr_group.store, mode="r"
        )
    except zarr.errors.GroupNotFoundError:
        logger.info("No variants in cache")
        cached_positions = pl.DataFrame({"chr_name": [], "chr_pos": []})
    else:
        logger.info(f"{meta_root} exists, getting cached positions")
        cached_positions = pl.DataFrame(
            {"chr_name": meta_root["chr_name"][:], "chr_pos": meta_root["chr_pos"][:]}
        )

    logger.info(f"Found {cached_positions.shape[0]} variants in cache")
    return cached_positions


class ZarrSampleMetadata(RootModel):
    root: list[str]

    def __len__(self) -> int:
        return len(self.root)
