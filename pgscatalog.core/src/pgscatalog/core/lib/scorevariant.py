"""This module contains classes that compose a ScoreVariant: a variant in a PGS
Catalog Scoring File."""

from enum import Enum
from functools import cached_property
from typing import Optional, ClassVar
from typing_extensions import Self
from pydantic import BaseModel, Field, model_validator, computed_field, field_validator


class Allele(BaseModel):
    """A class that represents an allele found in PGS Catalog scoring files
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
    """

    allele: str
    _valid_snp_bases: ClassVar[frozenset[str]] = frozenset({"A", "C", "T", "G"})

    @computed_field
    @property
    def is_snp(self) -> bool:
        """SNPs are the most common type of effect allele in PGS Catalog scoring
        files. More complex effect alleles, like HLAs or APOE genes, often require
        extra work to represent in genomes. Users should be warned about complex
        effect alleles.
        """
        return not frozenset(self.allele) - self._valid_snp_bases

    def __str__(self):
        return self.allele

    def __eq__(self, other):
        if isinstance(other, Allele):
            return self.allele == other.allele
        return False

    def __hash__(self):
        return hash(self.allele)


class EffectType(Enum):
    """Enums that represent inheritance models. The vast majority of variants have
    an additive effect type.

    This changes downstream PGS calculation:

    * ScoreVariants with an additive effect type will always be added to the PGS sum.
    * ScoreVariants with a dominant effect type are only added to the PGS sum if there is at least one copy of the effect allele.
    * ScoreVariants with a recessive effect type are only added to the PGS sum if there are two copies of the effect allele.

    >>> EffectType.ADDITIVE
    EffectType.ADDITIVE
    >>> str(EffectType.ADDITIVE)
    'additive'
    """

    RECESSIVE = "recessive"
    DOMINANT = "dominant"
    ADDITIVE = "additive"
    NONADDITIVE = "nonadditive"

    def __str__(self):
        return str(self.value)

    def __repr__(self):
        # pasting __repr__ output should be sufficient to construct the class
        return f"{type(self).__name__}.{self.name}"


class ScoreVariant(BaseModel):
    """A model representing a row from a PGS Catalog scoring file
    https://www.pgscatalog.org/downloads/#scoring_columns

    >>> variant = ScoreVariant(**{"rsID": None, "chr_name": "1", "chr_position": 1, "effect_allele": "A", "effect_weight": 0.5})
    >>> variant  # doctest: +ELLIPSIS
    ScoreVariant(rsID=None, chr_name='1', chr_position=1, effect_allele='A', ...
    >>> variant.is_complex
    False
    >>> variant.is_non_additive
    False
    >>> variant.is_harmonised
    False
    >>> variant.effect_type
    EffectType.ADDITIVE

    >>> variant_missing_positions = ScoreVariant(**{"rsID": None, "chr_name": None, "chr_position": None, "effect_allele": "A", "effect_weight": 0.5}) # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    TypeError: Bad position: rsid=None, chr_name=None, chr_position=None

    >>> harmonised_variant = ScoreVariant(**{"rsID": None, "chr_name": "1", "chr_position": 1, "effect_allele": "A", "effect_weight": 0.5, "hm_chr": "1", "hm_pos": 1, "hm_rsID": "rs1921", "hm_source": "ENSEMBL"})
    >>> harmonised_variant.is_harmonised
    True

    >>> variant_badly_harmonised = ScoreVariant(**{"rsID": None, "chr_name": "1", "chr_position": 1, "effect_allele": "A", "effect_weight": 0.5, "hm_chr": None, "hm_pos": None, "hm_rsID": "rs1921", "hm_source": "ENSEMBL"})
    Traceback (most recent call last):
    ...
    TypeError: Missing harmonised column data: hm_chr

    >>> variant_nonadditive = ScoreVariant(**{"rsID": None, "chr_name": "1", "chr_position": 1, "effect_allele": "A", "effect_weight": 0.5, "dosage_0_weight": 0, "dosage_1_weight": 1})
    >>> variant_nonadditive.is_non_additive
    True
    >>> variant_nonadditive.is_complex
    False
    >>> variant_nonadditive.effect_type
    EffectType.NONADDITIVE

    >>> variant_complex = ScoreVariant(**{"rsID": None, "chr_name": "1", "chr_position": 1, "effect_allele": "A", "effect_weight": 0.5, "is_haplotype": True})
    >>> variant_complex.is_complex
    True
    """

    # variant description
    rsID: Optional[str] = Field(
        default=None,
        title="dbSNP Accession ID (rsID)",
        description="The SNP’s rs ID. This column also contains HLA alleles in the standard notation (e.g. HLA-DQA1*0102) that aren’t always provided with chromosomal positions.",
    )
    chr_name: Optional[str] = Field(
        default=None,
        title="Location - Chromosome ",
        description="Chromosome name/number associated with the variant.",
    )
    chr_position: Optional[int] = Field(
        default=None,
        title="Location within the Chromosome",
        description="Chromosomal position associated with the variant.",
        gt=0,
    )
    effect_allele: Optional[str] = Field(
        default=None,
        title="Effect Allele",
        description="The allele that's dosage is counted (e.g. {0, 1, 2}) and multiplied by the variant's weight (effect_weight) when calculating score. The effect allele is also known as the 'risk allele'. Note: this does not necessarily need to correspond to the minor allele/alternative allele.",
    )
    other_allele: Optional[str] = Field(
        default=None,
        title="Other allele(s)",
        description="The other allele(s) at the loci. Note: this does not necessarily need to correspond to the reference allele.",
    )
    locus_name: Optional[str] = Field(
        default=None,
        title="Locus Name",
        description="This is kept in for loci where the variant may be referenced by the gene (APOE e4). It is also common (usually in smaller PGS) to see the variants named according to the genes they impact.",
    )
    is_haplotype: Optional[bool] = Field(
        default=False,
        title="FLAG: Haplotype",
        description="This is a TRUE/FALSE variable that flags whether the effect allele is a haplotype/diplotype rather than a single SNP. Constituent SNPs in the haplotype are semi-colon separated.",
    )
    is_diplotype: Optional[bool] = Field(
        default=False,
        title="FLAG: Diplotype",
        description="This is a TRUE/FALSE variable that flags whether the effect allele is a haplotype/diplotype rather than a single SNP. Constituent SNPs in the haplotype are semi-colon separated.",
    )
    imputation_method: Optional[str] = Field(
        default=None,
        title="Imputation Method",
        description="This described whether the variant was specifically called with a specific imputation or variant calling method. This is mostly kept to describe HLA-genotyping methods (e.g. flag SNP2HLA, HLA*IMP) that gives alleles that are not referenced by genomic position.",
    )
    variant_description: Optional[str] = Field(
        default=None,
        title="Variant Description",
        description="This field describes any extra information about the variant (e.g. how it is genotyped or scored) that cannot be captured by the other fields.",
    )
    inclusion_criteria: Optional[str] = Field(
        default=None,
        title="Score Inclusion Criteria",
        description="Explanation of when this variant gets included into the PGS (e.g. if it depends on the results from other variants).",
    )

    # weight information
    effect_weight: float = Field(
        title="Variant Weight",
        description="Value of the effect that is multiplied by the dosage of the effect allele (effect_allele) when calculating the score. Additional information on how the effect_weight was derived is in the weight_type field of the header, and score development method in the metadata downloads.",
    )
    is_interaction: Optional[bool] = Field(
        default=False,
        title="FLAG: Interaction",
        description="This is a TRUE/FALSE variable that flags whether the weight should be multiplied with the dosage of more than one variant. Interactions are demarcated with a _x_ between entries for each of the variants present in the interaction.",
    )
    is_dominant: Optional[bool] = Field(
        default=False,
        title="FLAG: Dominant Inheritance Model",
        description="This is a TRUE/FALSE variable that flags whether the weight should be added to the PGS sum if there is at least 1 copy of the effect allele (e.g. it is a dominant allele).",
    )
    is_recessive: Optional[bool] = Field(
        default=False,
        title="FLAG: Recessive Inheritance Model",
        description="This is a TRUE/FALSE variable that flags whether the weight should be added to the PGS sum only if there are 2 copies of the effect allele (e.g. it is a recessive allele).",
    )
    dosage_0_weight: Optional[float] = Field(
        default=None,
        title="Effect weight with 0 copy of the effect allele",
        description="Weights that are specific to different dosages of the effect_allele (e.g. {0, 1, 2} copies) can also be reported when the the contribution of the variants to the score is not encoded as additive, dominant, or recessive. In this case three columns are added corresponding to which variant weight should be applied for each dosage, where the column name is formated as dosage_#_weight where the # sign indicates the number of effect_allele copies.",
    )
    dosage_1_weight: Optional[float] = Field(
        default=None,
        title="Effect weight with 1 copy of the effect allele",
        description="Weights that are specific to different dosages of the effect_allele (e.g. {0, 1, 2} copies) can also be reported when the the contribution of the variants to the score is not encoded as additive, dominant, or recessive. In this case three columns are added corresponding to which variant weight should be applied for each dosage, where the column name is formated as dosage_#_weight where the # sign indicates the number of effect_allele copies.",
    )
    dosage_2_weight: Optional[float] = Field(
        default=None,
        title="Effect weight with 2 copies of the effect allele",
        description="Weights that are specific to different dosages of the effect_allele (e.g. {0, 1, 2} copies) can also be reported when the the contribution of the variants to the score is not encoded as additive, dominant, or recessive. In this case three columns are added corresponding to which variant weight should be applied for each dosage, where the column name is formated as dosage_#_weight where the # sign indicates the number of effect_allele copies.",
    )

    # other information
    OR: Optional[float] = Field(
        default=None,
        title="Odds Ratio",
        description="Author-reported effect sizes can be supplied to the Catalog. If no other effect_weight is given the weight is calculated using the log(OR) or log(HR).",
    )
    HR: Optional[float] = Field(
        default=None,
        title="Hazard Ratio",
        description="Author-reported effect sizes can be supplied to the Catalog. If no other effect_weight is given the weight is calculated using the log(OR) or log(HR).",
    )
    allelefrequency_effect: Optional[float] = Field(
        default=None,
        title="Effect Allele Frequency",
        description="Reported effect allele frequency, if the associated locus is a haplotype then haplotype frequency will be extracted.",
    )
    allelefrequency_effect_Ancestry: Optional[float] = Field(
        default=None,
        title="Population-specific effect allele frequency",
        description="Reported effect allele frequency in a specific population (described by the authors).",
    )

    # harmonised files - additional columns
    hm_source: Optional[str] = Field(
        default=None,
        title="Provider of the harmonized variant information",
        description="Data source of the variant position. Options include: ENSEMBL, liftover, author-reported (if being harmonized to the same build).",
    )
    hm_rsID: Optional[str] = Field(
        default=None,
        title="Harmonized rsID",
        description="Current rsID. Differences between this column and the author-reported column (rsID) indicate variant merges and annotation updates from dbSNP.",
    )
    hm_chr: Optional[str] = Field(
        default=None,
        title="Harmonized chromosome name",
        description="Chromosome that the harmonized variant is present on, preferring matches to chromosomes over patches present in later builds.",
    )
    hm_pos: Optional[int] = Field(
        default=None,
        title="Harmonized chromosome position",
        description="Chromosomal position (base pair location) where the variant is located, preferring matches to chromosomes over patches present in later builds.",
    )
    # this is a str on purpose
    hm_inferOtherAllele: Optional[str] = Field(
        default=None,
        title="Harmonized other alleles",
        description="If only the effect_allele is given we attempt to infer the non-effect/other allele(s) using Ensembl/dbSNP alleles.",
    )
    hm_match_chr: Optional[bool] = Field(
        default=None,
        title="FLAG: matching chromosome name",
        description="Used for QC. Only provided if the scoring file is being harmonized to the same genome build, and where the chromosome name is provided in the column chr_name.",
    )
    hm_match_pos: Optional[bool] = Field(
        default=None,
        title="FLAG: matching chromosome position",
        description="Used for QC. Only provided if the scoring file is being harmonized to the same genome build, and where the chromosome name is provided in the column chr_position.",
    )

    # helpful class attributes (not used by pydantic to instantiate a class)
    harmonised_columns: ClassVar[tuple[str]] = (
        "hm_rsID",
        "hm_chr",
        "hm_pos",
    )  # it's OK if (""hm_source", "hm_inferOtherAllele", "hm_match_chr", "hm_match_pos") are missing
    complex_columns: ClassVar[tuple[str]] = (
        "is_haplotype",
        "is_diplotype",
        "is_interaction",
    )
    non_additive_columns: ClassVar[tuple[str]] = (
        "dosage_0_weight",
        "dosage_1_weight",
        "dosage_2_weight",
    )

    # column names for output are used by __iter__ and when writing out
    output_fields: ClassVar[tuple[str]] = (
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

    @computed_field
    @cached_property
    def is_harmonised(self) -> bool:
        # simple check: do any of the harmonised columns have data?
        for x in self.harmonised_columns:
            if getattr(self, x) is not None:
                return True
        return False

    @computed_field
    @cached_property
    def is_complex(self) -> bool:
        # checking flag fields here, which are defaulted to False
        for x in self.complex_columns:
            if getattr(self, x):
                return True
        return False

    @computed_field
    @cached_property
    def is_non_additive(self) -> bool:
        # simple check: do any of the weight dosage columns have data?
        for x in self.non_additive_columns:
            if getattr(self, x) is not None:
                return True
        return False

    @computed_field
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

    @field_validator("effect_type")
    @classmethod
    def check_effect_type(cls, effect) -> EffectType:
        if not isinstance(effect, EffectType):
            raise TypeError("Bad effect type")
        return effect

    @model_validator(mode="after")
    def check_position(self) -> Self:
        rsid = getattr(self, "rsID", None)
        chr_name = getattr(self, "chr_name", None)
        chr_position = getattr(self, "chr_position", None)

        match (rsid, chr_name, chr_position):
            case None, str(), int():
                pass  # just genomic coordinates, good
            case str(), None, None:
                pass  # just rsID, good
            case str(), str(), int():
                pass  # all data are available, good
            case _:
                raise TypeError(f"Bad position: {rsid=}, {chr_name=}, {chr_position=} ")

        return self

    @model_validator(mode="after")
    def check_harmonised_columns(self) -> Self:
        if self.is_harmonised:
            for x in self.harmonised_columns:
                if getattr(self, x) is None:
                    raise TypeError(f"Missing harmonised column data: {x}")
        return self

    def __iter__(self):
        for attr in self.output_fields:
            yield getattr(self, attr)
