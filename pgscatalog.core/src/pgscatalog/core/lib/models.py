"""PGS Catalog pydantic models for data validation

Best way to reuse:

  * `from pgscatalog.core import models` and use `models.CatalogScoreVariant(**d)`

  * `import pgscatalog.core` and use fully qualified name: `pgscatalog.core.models.CatalogScoreVariant`)

"""

import enum
import itertools
import pathlib
from datetime import date
from functools import cached_property
from typing import ClassVar, Optional, Union, Any, Annotated

from typing_extensions import Self

from pydantic import (
    BaseModel,
    computed_field,
    model_serializer,
    Field,
    field_validator,
    model_validator,
    ConfigDict,
    field_serializer,
    RootModel,
)
from xopen import xopen

from pgscatalog.core.lib.genomebuild import GenomeBuild
from pgscatalog.core.lib.effecttype import EffectType


class Allele(BaseModel):
    """
    A class that represents an allele found in PGS Catalog scoring files

    >>> simple_ea = Allele(**{"allele": "A"})
    >>> simple_ea
    Allele(allele='A', is_snp=True)
    >>> str(simple_ea)
    'A'
    >>> Allele(**{"allele": "AG"})
    Allele(allele='AG', is_snp=True)
    >>> hla_example = Allele(**{"allele": "+"})
    >>> hla_example
    Allele(allele='+', is_snp=False)

    >>> Allele(allele="A")
    Allele(allele='A', is_snp=True)

    >>> Allele(allele="A/T").has_multiple_alleles
    True
    """

    allele: str
    _valid_snp_bases: ClassVar[frozenset[str]] = frozenset({"A", "C", "T", "G"})

    @computed_field  # type: ignore
    @cached_property
    def is_snp(self) -> bool:
        """SNPs are the most common type of effect allele in PGS Catalog scoring
        files. More complex effect alleles, like HLAs or APOE genes, often require
        extra work to represent in genomes. Users should be warned about complex
        effect alleles.
        """
        return not frozenset(self.allele) - self._valid_snp_bases

    @cached_property
    def has_multiple_alleles(self) -> bool:
        return "/" in self.allele

    @model_serializer(mode="plain", return_type=str)
    def serialize(self):
        """When dumping the model, flatten it to just return the allele as a string"""
        return self.allele

    def __str__(self):
        return self.allele

    def __eq__(self, other):
        if isinstance(other, Allele):
            return self.allele == other.allele
        return False

    def __hash__(self):
        return hash(self.allele)


class VariantType(str, enum.Enum):
    """Complex alleles are usually haplotypes/diplotypes and the gametic phase must be known to apply them accurately.

    See PGS Catalog Curation Guidelines: Appendix A – Special Cases for more information

    (Not supported by the calculator.)"""

    APOE_ALLELE = "APOE_allele"
    CYP_ALLELE = "CYP_allele"
    HLA_ALLELE = "HLA_allele"
    HLA_AA = "HLA_AA"
    HLA_SEROTYPE = "HLA_serotype"


class CatalogScoreVariant(BaseModel):
    """A model representing a row from a PGS Catalog scoring file, defined here:

    https://www.pgscatalog.org/downloads/#scoring_columns

    Implementation notes:

        - You should instantiate effect weight fields with strings (e.g. with csv.reader, which returns data as a list of strings)

        - The model always handles effect weights internally as strings and will coerce numeric input to strings when instantiated

        - Our string obsession comes from a desire to faithfully reproduce author submitted data and avoid introducing precision errors

    Extra / dynamically named fields:

    Only one type of dynamic field is supported. Ancestry specific allele frequency information uses labels defined by authors.

    An example from the first row from PGS000662:

    >>> variant_with_allelefrequency = {"chr_name": "1", "chr_position": 5743196, "effect_allele": "T", "other_allele": "C", "effect_weight": 0.102298257, "allelefrequency_effect_European": 0.067, "allelefrequency_effect_African": 0.439, "allelefrequency_effect_Asian": 0.113, "allelefrequency_effect_Hispanic": 0.157}
    >>> CatalogScoreVariant(**variant_with_allelefrequency)  # doctest: +ELLIPSIS
    CatalogScoreVariant(rsID=None, chr_name='1', chr_position=5743196..., allelefrequency_effect_European=0.067, allelefrequency_effect_African=0.439, allelefrequency_effect_Asian=0.113, allelefrequency_effect_Hispanic=0.157, ...)

    Extra field names which don't follow the pattern "allelefrequency_effect_{label}" will raise a ValueError:

    >>> bad_extra_fields = variant_with_allelefrequency | {"favourite_ice_cream": "vanilla"}
    >>> CatalogScoreVariant(**bad_extra_fields)
    Traceback (most recent call last):
    ...
    pydantic_core._pydantic_core.ValidationError: 1 validation error for CatalogScoreVariant
      Value error, Invalid extra fields detected: ['favourite_ice_cream'] ...

    Complex alleles are represented a little differently:

    >>> complex_allele = {"chr_name": 19, "effect_allele": "APOE_e2", "effect_weight": -0.5, "locus_name": "APOE", "is_haplotype": True, "variant_type": "APOE_allele", "variant_description": None}
    >>> CatalogScoreVariant(**complex_allele)
    CatalogScoreVariant(rsID=None, chr_name='19', chr_position=None, effect_allele=Allele(allele='APOE_e2', is_snp=False), other_allele=None, locus_name='APOE', is_haplotype=True, is_diplotype=False, imputation_method=None, variant_description=None, inclusion_criteria=None, effect_weight='-0.5', is_interaction=False, is_dominant=False, is_recessive=False, dosage_0_weight=None, dosage_1_weight=None, dosage_2_weight=None, OR=None, HR=None, allelefrequency_effect=None, hm_source=None, hm_rsID=None, hm_chr=None, hm_pos=None, hm_inferOtherAllele=None, hm_match_chr=None, hm_match_pos=None, variant_type=<VariantType.APOE_ALLELE: 'APOE_allele'>, variant_id='19::APOE_e2:', is_harmonised=False, is_complex=True, is_non_additive=False, effect_type=EffectType.ADDITIVE)

    Although effect weights are typed as optional, if all effect weight fields are missing then a model validator will raise a validation error:

    >>> CatalogScoreVariant(**{"chr_name": "19", "chr_position": 1, "effect_allele": "A", "effect_weight": None})  # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    pydantic_core._pydantic_core.ValidationError: 1 validation error for CatalogScoreVariant
      Value error, All effect weight fields are missing ...

    effect_weight can be missing if dosage_n_weight (non-additive) fields are all present

    However, dosage_n_weight fields must _all_ be present, if they're present:

    >>> CatalogScoreVariant(**{"chr_name": "19", "chr_position": 1, "effect_allele": "A", "dosage_0_weight": 0.1, "dosage_1_weight": None, "dosage_2_weight": None})  # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    pydantic_core._pydantic_core.ValidationError: 1 validation error for CatalogScoreVariant
      Value error, Dosage missing effect weight ...

    A variant may have all effect weight fields. During normalisation the standard effect_weight column will be used:

    >>> CatalogScoreVariant(**{"chr_name": "19", "chr_position": 1, "effect_allele": "A", "effect_weight": 0.05, "dosage_0_weight": 0, "dosage_1_weight": 0.1, "dosage_2_weight": 0.3})
    CatalogScoreVariant(rsID=None, chr_name='19', ..., is_non_additive=False, ...

    Note that is_non_additive is false if effect_weight column exists, although non-additive fields do exist.
    """

    model_config = ConfigDict(
        extra="allow"
    )  # extra fields are checked by a model validator

    # variant description
    rsID: Annotated[
        Optional[str],
        Field(
            default=None,
            title="dbSNP Accession ID (rsID)",
            description="The SNP’s rs ID. This column also contains HLA alleles in the standard notation (e.g. HLA-DQA1*0102) that aren’t always provided with chromosomal positions.",
        ),
    ]
    chr_name: Annotated[
        Optional[str],
        Field(
            default=None,
            title="Location - Chromosome ",
            description="Chromosome name/number associated with the variant.",
            coerce_numbers_to_str=True,
        ),
    ]
    chr_position: Annotated[
        Optional[int],
        Field(
            default=None,
            title="Location within the Chromosome",
            description="Chromosomal position associated with the variant.",
            gt=0,
        ),
    ]
    effect_allele: Annotated[
        Optional[Allele],
        Field(
            default=None,
            title="Effect Allele",
            description="The allele that's dosage is counted (e.g. {0, 1, 2}) and multiplied by the variant's weight (effect_weight) when calculating score. The effect allele is also known as the 'risk allele'. Note: this does not necessarily need to correspond to the minor allele/alternative allele.",
        ),
    ]
    other_allele: Annotated[
        Optional[Allele],
        Field(
            default=None,
            title="Other allele(s)",
            description="The other allele(s) at the loci. Note: this does not necessarily need to correspond to the reference allele.",
        ),
    ]
    locus_name: Annotated[
        Optional[str],
        Field(
            default=None,
            title="Locus Name",
            description="This is kept in for loci where the variant may be referenced by the gene (APOE e4). It is also common (usually in smaller PGS) to see the variants named according to the genes they impact.",
        ),
    ]
    is_haplotype: Annotated[
        Optional[bool],
        Field(
            default=False,
            title="FLAG: Haplotype",
            description="This is a TRUE/FALSE variable that flags whether the effect allele is a haplotype/diplotype rather than a single SNP. Constituent SNPs in the haplotype are semi-colon separated.",
        ),
    ]
    is_diplotype: Annotated[
        Optional[bool],
        Field(
            default=False,
            title="FLAG: Diplotype",
            description="This is a TRUE/FALSE variable that flags whether the effect allele is a haplotype/diplotype rather than a single SNP. Constituent SNPs in the haplotype are semi-colon separated.",
        ),
    ]
    imputation_method: Annotated[
        Optional[str],
        Field(
            default=None,
            title="Imputation Method",
            description="This described whether the variant was specifically called with a specific imputation or variant calling method. This is mostly kept to describe HLA-genotyping methods (e.g. flag SNP2HLA, HLA*IMP) that gives alleles that are not referenced by genomic position.",
        ),
    ]
    variant_description: Annotated[
        Optional[str],
        Field(
            default=None,
            title="Variant Description",
            description="This field describes any extra information about the variant (e.g. how it is genotyped or scored) that cannot be captured by the other fields.",
        ),
    ]
    inclusion_criteria: Annotated[
        Optional[str],
        Field(
            default=None,
            title="Score Inclusion Criteria",
            description="Explanation of when this variant gets included into the PGS (e.g. if it depends on the results from other variants).",
        ),
    ]

    # weight information
    # all effect weight fields are handled as strings internally to faithfully reproduce author-uploaded scores (i.e. avoid any floating point errors)
    # weight fields are validated by effect_weight_must_float and check_effect_weights
    # note: string coercion only happens if class is instantiated with numeric data
    # to avoid precision errors, always instantiate with strings, e.g. csv.reader which always return a list of strings
    effect_weight: Annotated[
        Optional[str],
        Field(
            default=None,
            title="Variant Weight",
            description="Value of the effect that is multiplied by the dosage of the effect allele (effect_allele) when calculating the score. Additional information on how the effect_weight was derived is in the weight_type field of the header, and score development method in the metadata downloads.",
            coerce_numbers_to_str=True,
        ),
    ]
    is_interaction: Annotated[
        Optional[bool],
        Field(
            default=False,
            title="FLAG: Interaction",
            description="This is a TRUE/FALSE variable that flags whether the weight should be multiplied with the dosage of more than one variant. Interactions are demarcated with a _x_ between entries for each of the variants present in the interaction.",
        ),
    ]
    is_dominant: Annotated[
        Optional[bool],
        Field(
            default=False,
            title="FLAG: Dominant Inheritance Model",
            description="This is a TRUE/FALSE variable that flags whether the weight should be added to the PGS sum if there is at least 1 copy of the effect allele (e.g. it is a dominant allele).",
        ),
    ]
    is_recessive: Annotated[
        Optional[bool],
        Field(
            default=False,
            title="FLAG: Recessive Inheritance Model",
            description="This is a TRUE/FALSE variable that flags whether the weight should be added to the PGS sum only if there are 2 copies of the effect allele (e.g. it is a recessive allele).",
        ),
    ]
    dosage_0_weight: Annotated[
        Optional[str],
        Field(
            default=None,
            title="Effect weight with 0 copy of the effect allele",
            description="Weights that are specific to different dosages of the effect_allele (e.g. {0, 1, 2} copies) can also be reported when the the contribution of the variants to the score is not encoded as additive, dominant, or recessive. In this case three columns are added corresponding to which variant weight should be applied for each dosage, where the column name is formated as dosage_#_weight where the # sign indicates the number of effect_allele copies.",
            coerce_numbers_to_str=True,
        ),
    ]
    dosage_1_weight: Annotated[
        Optional[str],
        Field(
            default=None,
            title="Effect weight with 1 copy of the effect allele",
            description="Weights that are specific to different dosages of the effect_allele (e.g. {0, 1, 2} copies) can also be reported when the the contribution of the variants to the score is not encoded as additive, dominant, or recessive. In this case three columns are added corresponding to which variant weight should be applied for each dosage, where the column name is formated as dosage_#_weight where the # sign indicates the number of effect_allele copies.",
            coerce_numbers_to_str=True,
        ),
    ]
    dosage_2_weight: Annotated[
        Optional[str],
        Field(
            default=None,
            title="Effect weight with 2 copies of the effect allele",
            description="Weights that are specific to different dosages of the effect_allele (e.g. {0, 1, 2} copies) can also be reported when the the contribution of the variants to the score is not encoded as additive, dominant, or recessive. In this case three columns are added corresponding to which variant weight should be applied for each dosage, where the column name is formated as dosage_#_weight where the # sign indicates the number of effect_allele copies.",
            coerce_numbers_to_str=True,
        ),
    ]

    # other information
    OR: Annotated[
        Optional[float],
        Field(
            default=None,
            title="Odds Ratio",
            description="Author-reported effect sizes can be supplied to the Catalog. If no other effect_weight is given the weight is calculated using the log(OR) or log(HR).",
        ),
    ]
    HR: Annotated[
        Optional[float],
        Field(
            default=None,
            title="Hazard Ratio",
            description="Author-reported effect sizes can be supplied to the Catalog. If no other effect_weight is given the weight is calculated using the log(OR) or log(HR).",
        ),
    ]
    allelefrequency_effect: Annotated[
        Optional[float],
        Field(
            default=None,
            title="Effect Allele Frequency",
            description="Reported effect allele frequency, if the associated locus is a haplotype then haplotype frequency will be extracted.",
            ge=0,
        ),
    ]

    # harmonised files - additional columns
    hm_source: Annotated[
        Optional[str],
        Field(
            default=None,
            title="Provider of the harmonized variant information",
            description="Data source of the variant position. Options include: ENSEMBL, liftover, author-reported (if being harmonized to the same build).",
        ),
    ]
    hm_rsID: Annotated[
        Optional[str],
        Field(
            default=None,
            title="Harmonized rsID",
            description="Current rsID. Differences between this column and the author-reported column (rsID) indicate variant merges and annotation updates from dbSNP.",
        ),
    ]
    hm_chr: Annotated[
        Optional[str],
        Field(
            default=None,
            title="Harmonized chromosome name",
            description="Chromosome that the harmonized variant is present on, preferring matches to chromosomes over patches present in later builds.",
        ),
    ]
    hm_pos: Annotated[
        Optional[int],
        Field(
            ge=0,
            default=None,
            title="Harmonized chromosome position",
            description="Chromosomal position (base pair location) where the variant is located, preferring matches to chromosomes over patches present in later builds.",
        ),
    ]
    hm_inferOtherAllele: Annotated[
        Optional[Allele],
        Field(
            default=None,
            title="Harmonized other alleles",
            description="If only the effect_allele is given we attempt to infer the non-effect/other allele(s) using Ensembl/dbSNP alleles.",
        ),
    ]
    hm_match_chr: Annotated[
        Optional[bool],
        Field(
            default=None,
            title="FLAG: matching chromosome name",
            description="Used for QC. Only provided if the scoring file is being harmonized to the same genome build, and where the chromosome name is provided in the column chr_name.",
        ),
    ]
    hm_match_pos: Annotated[
        Optional[bool],
        Field(
            default=None,
            title="FLAG: matching chromosome position",
            description="Used for QC. Only provided if the scoring file is being harmonized to the same genome build, and where the chromosome name is provided in the column chr_position.",
        ),
    ]
    variant_type: Annotated[
        Optional[VariantType],
        Field(
            default=None,
            title="Complex alleles only: how is the variant name formatted?",
        ),
    ]

    # helpful class attributes (not used by pydantic to instantiate a class)
    harmonised_columns: ClassVar[tuple[str, str, str]] = ("hm_rsID", "hm_chr", "hm_pos")
    complex_columns: ClassVar[tuple[str, str, str]] = (
        "is_haplotype",
        "is_diplotype",
        "is_interaction",
    )
    non_additive_columns: ClassVar[tuple[str, str, str]] = (
        "dosage_0_weight",
        "dosage_1_weight",
        "dosage_2_weight",
    )

    # these computed fields provide convenient properties e.g. when normalising variants
    @computed_field  # type: ignore
    @cached_property
    def variant_id(self) -> str:
        """ID = chr:pos:effect_allele:other_allele"""
        return ":".join(
            [
                str(getattr(self, k) or "")  # correctly handles None elements
                for k in ["chr_name", "chr_position", "effect_allele", "other_allele"]
            ]
        )

    @computed_field  # type: ignore
    @cached_property
    def is_harmonised(self) -> bool:
        # simple check: do any of the harmonised columns have data?
        for x in self.harmonised_columns:
            if getattr(self, x) is not None:
                return True
        return False

    @computed_field  # type: ignore
    @cached_property
    def is_complex(self) -> bool:
        is_complex = not getattr(self.effect_allele, "is_snp", False)
        for x in self.complex_columns:
            if getattr(self, x):
                is_complex = True

        return is_complex

    @computed_field  # type: ignore
    @cached_property
    def is_non_additive(self) -> bool:
        if self.effect_weight is not None:
            # if there's an effect weight value, we can work with it
            non_additive = False
        else:
            # dosage columns are trickier
            non_additive = any(
                getattr(self, col) is not None for col in self.non_additive_columns
            )
        return non_additive

    @computed_field  # type: ignore
    @cached_property
    def effect_type(self) -> EffectType:
        match (self.is_recessive, self.is_dominant, self.is_non_additive):
            case False, False, False:
                effect = EffectType.ADDITIVE
            case True, False, False:
                effect = EffectType.RECESSIVE
            case False, True, False:
                effect = EffectType.DOMINANT
            case False, False, True:
                effect = EffectType.NONADDITIVE
            case _:
                raise ValueError(
                    f"Can't determine effect type: {self.is_recessive=}, {self.is_dominant=}, {self.is_non_additive=}"
                )

        return effect

    # start validating individual fields from here (happens before model validation)
    # ordered by mode: "before" field validators work with raw data (e.g. strings/ints)
    # "after" field validators work with fields that are the correct type
    @field_validator(
        "effect_allele", "other_allele", "hm_inferOtherAllele", mode="before"
    )
    @classmethod
    def alleles_must_parse(cls, value: Any) -> Allele:
        if isinstance(value, str):
            return Allele(allele=value)
        else:
            raise ValueError(f"Can't parse {value=}")

    @field_validator(
        "rsID",
        "chr_name",
        "chr_position",
        "is_haplotype",
        "is_diplotype",
        "hm_chr",
        "hm_pos",
        "hm_match_chr",
        "hm_match_pos",
        "allelefrequency_effect",
        "variant_type",
        mode="before",
    )
    @classmethod
    def empty_string_to_none(cls, v: Any) -> Optional[Any]:
        if isinstance(v, str) and v.strip() == "":
            return None
        return v

    @field_validator("rsID", mode="after")
    @classmethod
    def check_rsid_format(cls, rsid: Optional[str]) -> Optional[str]:
        if rsid == ".":
            rsid = None

        if (
            rsid is None
            or rsid.startswith("rs")
            or rsid.startswith("ss")
            or rsid.startswith("HLA")
        ):
            return rsid
        else:
            raise ValueError("rsid field must start with rs or ss or HLA")

    @field_validator(
        "effect_weight",
        "dosage_0_weight",
        "dosage_1_weight",
        "dosage_2_weight",
        mode="after",
    )
    @classmethod
    def effect_weight_must_float(cls, weight: Optional[str]) -> Optional[str]:
        # reminder: weight fields default to None because check_effect_weights validates the model after instantiation
        # default values for fields are never validated
        # but any non-default value is validated, and must be coercible to float
        if weight is not None:
            _ = float(weight)  # will raise a ValueError if conversion fails

        return weight

    # start validating the entire model now (i.e. after instantiation, all fields are set)
    @model_validator(mode="after")
    def check_extra_fields(self) -> Self:
        """Only allelefrequency_effect_{ancestry} is supported as an extra field
        {ancestry} is dynamic and set by submitters"""
        if self.model_extra is not None:
            extra: list[str] = list(self.model_extra.keys())

            if extra:
                field_match: list[bool] = [
                    x.startswith("allelefrequency_effect_") for x in extra
                ]
                if not all(field_match):
                    bad_extra_fields: list[str] = list(
                        itertools.compress(extra, [not x for x in field_match])
                    )
                    raise ValueError(
                        f"Invalid extra fields detected: {bad_extra_fields}"
                    )
                else:
                    for field in extra:
                        # make sure allele frequency is a float or raise a value error
                        allelefrequency: float = float(getattr(self, field))
                        setattr(self, field, allelefrequency)

        return self

    @model_validator(mode="after")
    def check_effect_weights(self) -> Self:
        match (
            self.effect_weight,
            self.dosage_0_weight,
            self.dosage_1_weight,
            self.dosage_2_weight,
        ):
            case None, None, None, None:
                raise ValueError("All effect weight fields are missing")
            case None, zero, one, two if any(x is None for x in (zero, one, two)):
                raise ValueError("Dosage missing effect weight")
            case _:
                return self

    @model_validator(mode="after")
    def check_position(self) -> Self:
        match (self.rsID, self.chr_name, self.chr_position):
            case str() | None, str(), int():
                # mandatory coordinates with optional rsid
                pass
            case str(), str() | None, str() | None:
                # mandatory rsid with optional coordinates
                pass
            case _:
                if self.is_complex:
                    # complex variants can have odd non-standard positions
                    # e.g. e2 allele of APOE gene
                    pass
                else:
                    raise ValueError(
                        f"Bad position: {self.rsID=}, {self.chr_name=}, {self.chr_position=} "
                        f"for variant {self=}"
                    )

        return self

    @model_validator(mode="after")
    def check_complex_variants(self) -> Self:
        if self.variant_type is not None:
            if self.is_haplotype is None and self.is_diplotype is None:
                raise ValueError(
                    f"Complex alleles must have is_haplotype or is_diplotype fields: {self}"
                )

            if getattr(self.effect_allele, "is_snp", False):
                raise ValueError(
                    f"Complex variants can't have SNP effect alleles: {self}"
                )

        return self


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
    pydantic_core._pydantic_core.ValidationError: 1 validation error for ScoreVariant
      Value error, Bad position: self.rsID=None, self.chr_name=None, self.chr_position=None...
      ...

    >>> harmonised_variant = ScoreVariant(**{"rsID": None, "chr_name": "1", "chr_position": 1, "effect_allele": "A", "effect_weight": 0.5, "hm_chr": "1", "hm_pos": 1, "hm_rsID": "rs1921", "hm_source": "ENSEMBL",  "row_nr": 0, "accession": "test"})
    >>> harmonised_variant.is_harmonised
    True

    >>> variant_nonadditive = ScoreVariant(**{"rsID": None, "chr_name": "1", "chr_position": 1, "effect_allele": "A", "dosage_0_weight": 0, "dosage_1_weight": 1,  "dosage_2_weight": 0, "row_nr": 0, "accession": "test"})
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

    model_config = ConfigDict(use_enum_values=True)

    row_nr: Annotated[
        int,
        Field(
            title="Row number",
            description="Row number of variant in scoring file (first variant = 0)",
            ge=0,
        ),
    ]
    accession: Annotated[
        str, Field(title="Accession", description="Accession of score variant")
    ]
    is_duplicated: Annotated[
        Optional[bool],
        Field(
            default=False,
            title="Duplicated variant",
            description="In a list of variants with the same accession, is ID duplicated?",
        ),
    ]

    # column names for output are used by __iter__ and when writing out
    output_fields: ClassVar[tuple[str, ...]] = (
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


class ScoreFormatVersion(str, enum.Enum):
    """See https://www.pgscatalog.org/downloads/#scoring_changes
    v1 was deprecated in December 2021
    """

    v2 = "2.0"


class ScoreHeader(BaseModel):
    """Headers store useful metadata about a scoring file.

    Data validation is less strict than the CatalogScoreHeader, to make
    it easier for people to use custom scoring files with the PGS Catalog Calculator.

    >>> ScoreHeader(**{"pgs_id": "PGS123456", "trait_reported": "testtrait", "genome_build": "GRCh38"})
    ScoreHeader(pgs_id='PGS123456', pgs_name=None, trait_reported='testtrait', genome_build=GenomeBuild.GRCh38)

    >>> from ._config import Config
    >>> testpath = Config.ROOT_DIR / "tests" / "data" / "PGS000001_hmPOS_GRCh38.txt.gz"
    >>> ScoreHeader.from_path(testpath).row_count
    77
    """

    pgs_id: Annotated[str, Field(title="PGS identifier")]
    pgs_name: Annotated[Optional[str], Field(description="PGS name", default=None)]
    trait_reported: Annotated[str, Field(description="Trait name")]
    # genome build is Optional because "NR" is represented internally as None
    genome_build: Annotated[Optional[GenomeBuild], Field(description="Genome build")]

    _path: Optional[pathlib.Path]

    @property
    def is_harmonised(self):
        # custom scores can't be harmonised o_o
        return False

    @field_validator("genome_build", mode="before")
    @classmethod
    def parse_genome_build(cls, value: str) -> Optional[GenomeBuild]:
        if value == "NR":
            return None
        else:
            return GenomeBuild.from_string(value)

    @field_serializer("genome_build")
    def serialize_genomebuild(self, genome_build, _info):
        return genome_build.value if genome_build is not None else "NR"

    @cached_property
    def row_count(self) -> int:
        """Calculate the number of variants in the scoring file by counting the number of rows"""
        path: Optional[pathlib.Path] = getattr(self, "_path", None)
        if path is None:
            raise TypeError("Can't calculate row count without path")
        else:
            with xopen(path, "rt") as fh:
                # skip header line with - 1
                n_variants = sum(1 for x in fh if not x.startswith("#")) - 1
                if n_variants == 0:
                    raise ValueError(f"No variants in {path}")

        return n_variants

    @classmethod
    def from_path(cls, path):
        header = {}

        def generate_header(f):
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#"):
                    if "=" in line:
                        yield line.strip()
                else:
                    # stop reading lines
                    break

        with xopen(path, "rt") as f:
            header_text = generate_header(f)

            for item in header_text:
                key, value = item.split("=")
                header[key[1:]] = value  # drop # character from key

        scoreheader = cls(**header)
        scoreheader._path = pathlib.Path(path)
        return scoreheader


class CatalogScoreHeader(ScoreHeader):
    """A ScoreHeader that validates the PGS Catalog Scoring File header standard

    https://www.pgscatalog.org/downloads/#dl_ftp_scoring

    >>> from ._config import Config
    >>> testpath = Config.ROOT_DIR / "tests" / "data" / "PGS000001_hmPOS_GRCh38.txt.gz"
    >>> test = CatalogScoreHeader.from_path(testpath) # doctest: +ELLIPSIS
    >>> test # doctest: +ELLIPSIS
    CatalogScoreHeader(pgs_id='PGS000001', pgs_name='PRS77_BC', trait_reported='Breast cancer', genome_build=None, format_version=<ScoreFormatVersion.v2: '2.0'>, trait_mapped=['breast carcinoma'], trait_efo=['EFO_0000305'], variants_number=77, weight_type=None, pgp_id='PGP000001', citation='Mavaddat N et al. J Natl Cancer Inst (2015). doi:10.1093/jnci/djv036', HmPOS_build=GenomeBuild.GRCh38, HmPOS_date=datetime.date(2022, 7, 29), HmPOS_match_pos='{"True": null, "False": null}', HmPOS_match_chr='{"True": null, "False": null}')
    >>> test.variants_number == test.row_count
    True
    """

    format_version: ScoreFormatVersion
    trait_mapped: Annotated[list[str], Field(description="Trait name")]
    trait_efo: Annotated[
        list[str], Field(description="Ontology trait name, e.g. 'breast carcinoma")
    ]
    variants_number: Annotated[
        int,
        Field(gt=0, description="Number of variants listed in the PGS", default=None),
    ]
    # note: we'll make sure to serialise None values here and in genome_build as string "NR"
    weight_type: Annotated[
        Optional[str], Field(description="Variant weight type", default=None)
    ]

    ##SOURCE INFORMATION
    pgp_id: str
    citation: str
    ##HARMONIZATION DETAILS
    HmPOS_build: Annotated[Optional[GenomeBuild], Field(default=None)]
    HmPOS_date: Annotated[Optional[date], Field(default=None)]
    HmPOS_match_pos: Annotated[Optional[str], Field(default=None)]
    HmPOS_match_chr: Annotated[Optional[str], Field(default=None)]

    # note: only included when different from default
    license: Annotated[
        Optional[str],
        Field(
            "PGS obtained from the Catalog should be cited appropriately, and "
            "used in accordance with any licensing restrictions set by the authors. See "
            "EBI Terms of Use (https://www.ebi.ac.uk/about/terms-of-use/) for additional "
            "details.",
            repr=False,
        ),
    ]

    @field_validator("trait_mapped", "trait_efo", mode="before")
    @classmethod
    def split_traits(cls, trait: str) -> list[str]:
        if isinstance(trait, str):
            traits = trait.split("|")
            if len(traits) == 0:
                raise ValueError("No traits defined")
            return traits
        raise ValueError(f"Can't parse trait string: {trait}")

    @classmethod
    def _check_accession(cls, value: str, prefix: str) -> str:
        if not value.startswith(prefix):
            raise ValueError(f"{value} doesn't start with {prefix}")
        if len(value) != 9:
            raise ValueError(f"Invalid accession format: {value}")
        return value

    @field_validator("pgs_id")
    @classmethod
    def check_pgs_id(cls, pgs_id: str) -> str:
        return cls._check_accession(pgs_id, "PGS")

    @field_validator("pgp_id")
    @classmethod
    def check_pgp_id(cls, pgp_id: str) -> str:
        return cls._check_accession(pgp_id, "PGP")

    @field_validator("genome_build", "HmPOS_build", mode="before")
    @classmethod
    def parse_genome_build(cls, value: str) -> Optional[GenomeBuild]:
        if value == "NR":
            return None
        else:
            return GenomeBuild.from_string(value)

    @field_serializer("genome_build", "HmPOS_build")
    def serialize_genomebuild(self, genome_build, _info):
        return genome_build.value if genome_build is not None else "NR"

    @field_validator("format_version")
    @classmethod
    def check_format_version(cls, version: ScoreFormatVersion) -> ScoreFormatVersion:
        if version != ScoreFormatVersion.v2:
            raise ValueError(f"Invalid format_version: {version}")
        return version

    @field_validator("weight_type", mode="after")
    @classmethod
    def parse_weight_type(cls, value: Optional[str]) -> Optional[str]:
        if value == "NR":
            value = None
        return value

    @property
    def is_harmonised(self):
        if self.HmPOS_build is None and self.HmPOS_date is None:
            return False
        else:
            return True


class VariantLog(BaseModel):
    """This model consists of variant-level statistics we need to summarise in the ScoreLog

    Can't just reuse ScoreVariants because failed harmonisation can create invalid ScoreVariants (e.g. missing genomic coordinates)

    If ScoreLogs are composed of ScoreVariants then the data would be revalidated on instantiation and raise ValidationErrors

    Instead, just create VariantLogs from a subset of ScoreVariant fields
    """

    hm_source: Optional[str] = None
    is_complex: bool


class ScoreLog(BaseModel):
    """A log that includes header information and variant summary statistics

    >>> header = CatalogScoreHeader(pgs_id='PGS000001', pgs_name='PRS77_BC', trait_reported='Breast cancer', genome_build=None, format_version=ScoreFormatVersion.v2, trait_mapped='breast carcinoma', trait_efo='EFO_0000305', variants_number=77, weight_type="NR", pgp_id='PGP000001', citation='Mavaddat N et al. J Natl Cancer Inst (2015). doi:10.1093/jnci/djv036', HmPOS_build="GRCh38", HmPOS_date="2022-07-29")
    >>> harmonised_variant = ScoreVariant(**{"rsID": None, "chr_name": "1", "chr_position": 1, "effect_allele": "HLA-DQ", "effect_weight": 0.5, "hm_chr": "1", "hm_pos": 1, "hm_rsID": "rs1921", "hm_source": "ENSEMBL",  "row_nr": 0, "accession": "test"})
    >>> variant_log = harmonised_variant.model_dump(include={"hm_source", "is_complex"})
    >>> scorelog = ScoreLog(header=header, compatible_effect_type=True, variant_logs=[VariantLog(**variant_log)])  # doctest: +ELLIPSIS
    >>> scorelog
    ScoreLog(header=CatalogScoreHeader(...), compatible_effect_type=True, has_complex_alleles=True, pgs_id='PGS000001', is_harmonised=True, sources=['ENSEMBL'])

    In the original scoring file header there were 77 variants:

    >>> scorelog.header.variants_number
    77

    But we've only got 1 ScoreVariant:

    >>> scorelog.n_actual_variants
    1
    >>> scorelog.variant_count_difference
    76
    >>> scorelog.variants_are_missing
    True

    Maybe they were all filtered out during normalisationIt's important to log and warn when this happens.

    >>> scorelog.sources
    ['ENSEMBL']

    >>> scorelog.model_dump()  # doctest: +ELLIPSIS
    {'header': {'pgs_id': 'PGS000001', ...}, 'compatible_effect_type': True, 'has_complex_alleles': True, 'pgs_id': 'PGS000001', 'is_harmonised': True, 'sources': ['ENSEMBL']}
    """

    model_config = ConfigDict(use_enum_values=True)

    header: Union[ScoreHeader, CatalogScoreHeader] = Field(
        description="Metadata from the scoring file header"
    )
    variant_logs: Optional[list[VariantLog]] = Field(
        description="A list of variants associated with the header. Some may be filtered out during normalisation.",
        exclude=True,
        repr=False,
    )
    compatible_effect_type: bool = Field(
        description="Did all variants in this score contain compatible effect types? (i.e. additive / recessive / dominant)"
    )

    @computed_field  # type: ignore
    @cached_property
    def has_complex_alleles(self) -> bool:
        """Do any variants contain complex alleles? e.g. HLA/APOE"""
        has_complex = False
        if self.variant_logs is not None:
            has_complex = any(x.is_complex for x in self.variant_logs)
        return has_complex

    @computed_field  # type: ignore
    @property
    def pgs_id(self) -> str:
        return self.header.pgs_id

    @computed_field  # type: ignore
    @property
    def is_harmonised(self) -> bool:
        return self.header.is_harmonised

    @computed_field  # type: ignore
    @cached_property
    def sources(self) -> Optional[list[str]]:
        unique_sources: Optional[list[str]] = None
        if self.variant_logs is not None:
            sources: list[Optional[str]] = [x.hm_source for x in self.variant_logs]
            filtered_sources = [x for x in sources if x is not None]
            if not filtered_sources:
                unique_sources = None
            else:
                unique_sources = list(set(filtered_sources))
        return unique_sources

    @property
    def n_actual_variants(self) -> Optional[int]:
        # this distinction is useful if variants have been filtered out
        if self.variant_logs is not None:
            return len(self.variant_logs)
        else:
            return None

    @cached_property
    def variant_count_difference(self) -> Optional[int]:
        # grab directly from header
        header_variants = getattr(self.header, "variants_number", None)
        if header_variants is None:
            # parse by counting the number of rows in the scoring file
            header_variants = self.header.row_count

        if self.n_actual_variants is None:
            variant_difference: Optional[int] = None
        else:
            variant_difference = abs(header_variants - self.n_actual_variants)

        return variant_difference

    @property
    def variants_are_missing(self) -> bool:
        return self.variant_count_difference != 0


class ScoreLogs(RootModel):
    """A container of ScoreLog to simplify serialising to a JSON list"""

    root: list[ScoreLog]
