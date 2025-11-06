"""
Pydantic models for loading and saving data from zarr attributes (metadata)

These models help to read and write structured data about variants and samples
This data is helpful when working with the genotype array (variants = row names,
samples = column names).
"""

from __future__ import annotations

import logging
from copy import deepcopy
from typing import Annotated

import polars as pl
from pydantic import AfterValidator, BaseModel, PositiveInt, model_validator

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
    sampleset: list[str]
    filename: list[str]
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


def is_valid_allele(alleles: list[list[str] | None] | list[str | None]):
    valid_alleles = {"A", "C", "T", "G"}
    for allele in alleles:
        if allele is not None:
            if not set(allele).issubset(valid_alleles):
                raise ValueError(f"Invalid allele: {allele}")
    return alleles


class ZarrSampleMetadata(BaseModel):
    samples: list[str]

    def __len__(self) -> int:
        return len(self.samples)
