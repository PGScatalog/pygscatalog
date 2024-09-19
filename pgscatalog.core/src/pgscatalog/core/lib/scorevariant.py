from typing import Optional, ClassVar, Literal
from pydantic import (
    Field,
    field_serializer,
)

from .effecttype import EffectType
from .models import CatalogScoreVariant


class ScoreVariant(CatalogScoreVariant):
    """This model includes attributes useful for processing and normalising variants

    >>> variant = ScoreVariant(**{"rsID": None, "chr_name": "1", "chr_position": 1, "effect_allele": "A", "effect_weight": 0.5, "row_nr": 0, "accession": "test"})
    >>> variant  # doctest: +ELLIPSIS
    ScoreVariant(rsID=None, chr_name='1', chr_position=1, effect_allele=Allele(allele='A', ...
    >>> variant.is_complex
    False
    >>> variant.is_non_additive
    False
    >>> variant.is_harmonised
    False
    >>> variant.effect_type
    EffectType.ADDITIVE

    >>> variant_missing_positions = ScoreVariant(**{"rsID": None, "chr_name": None, "chr_position": None, "effect_allele": "A", "effect_weight": 0.5,  "row_nr": 0, "accession": "test"}) # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    TypeError: Bad position: self.rsID=None, self.chr_name=None, self.chr_position=None

    >>> harmonised_variant = ScoreVariant(**{"rsID": None, "chr_name": "1", "chr_position": 1, "effect_allele": "A", "effect_weight": 0.5, "hm_chr": "1", "hm_pos": 1, "hm_rsID": "rs1921", "hm_source": "ENSEMBL",  "row_nr": 0, "accession": "test"})
    >>> harmonised_variant.is_harmonised
    True

    >>> variant_badly_harmonised = ScoreVariant(**{"rsID": None, "chr_name": "1", "chr_position": 1, "effect_allele": "A", "effect_weight": 0.5, "hm_chr": None, "hm_pos": None, "hm_rsID": "rs1921", "hm_source": "ENSEMBL",  "row_nr": 0, "accession": "test"})
    Traceback (most recent call last):
    ...
    TypeError: Missing harmonised column data: hm_chr

    >>> variant_nonadditive = ScoreVariant(**{"rsID": None, "chr_name": "1", "chr_position": 1, "effect_allele": "A", "effect_weight": 0.5, "dosage_0_weight": 0, "dosage_1_weight": 1,  "row_nr": 0, "accession": "test"})
    >>> variant_nonadditive.is_non_additive
    True
    >>> variant_nonadditive.is_complex
    False
    >>> variant_nonadditive.effect_type
    EffectType.NONADDITIVE

    >>> variant_complex = ScoreVariant(**{"rsID": None, "chr_name": "1", "chr_position": 1, "effect_allele": "A", "effect_weight": 0.5, "is_haplotype": True,  "row_nr": 0, "accession": "test"})
    >>> variant_complex.is_complex
    True
    """

    row_nr: int = Field(
        title="Row number",
        description="Row number of variant in scoring file (first variant = 0)",
    )
    accession: str = Field(title="Accession", description="Accession of score variant")
    is_duplicated: Optional[bool] = Field(
        default=False,
        title="Duplicated variant",
        description="In a list of variants with the same accession, is ID duplicated?",
    )

    # column names for output are used by __iter__ and when writing out
    output_fields: ClassVar[
        tuple[
            Literal["chr_name"],
            Literal["chr_position"],
            Literal["effect_allele"],
            Literal["other_allele"],
            Literal["effect_weight"],
            Literal["effect_type"],
            Literal["is_duplicated"],
            Literal["accession"],
            Literal["row_nr"],
        ]
    ] = (
        "chr_name",
        "chr_position",
        "effect_allele",
        "other_allele",
        "effect_weight",
        "effect_type",
        "is_duplicated",
        "accession",
        "row_nr",
    )

    def __iter__(self):
        for attr in self.output_fields:
            yield getattr(self, attr)

    @field_serializer("effect_type")
    def serialize_effect_type(self, effect_type: EffectType) -> str:
        """Convert enum to string during serialisation"""
        return effect_type.value
