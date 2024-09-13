""" PGS Catalog pydantic models for data validation

Best way to reuse:

  * `from pgscatalog.core import models` and use `models.CatalogScoreVariant(**d)`

  * `import pgscatalog.core` and use fully qualified name: `pgscatalog.core.models.CatalogScoreVariant`)

"""
import enum
from datetime import date
from functools import cached_property
from typing import ClassVar, Optional
from typing_extensions import Self

from pydantic import (
    BaseModel,
    computed_field,
    model_serializer,
    Field,
    field_validator,
    model_validator,
)
from xopen import xopen

from ..lib import EffectType, GenomeBuild


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

    >>> Allele(allele="A")
    Allele(allele='A', is_snp=True)

    >>> Allele(allele="A/T").has_multiple_alleles
    True
    """

    allele: str
    _valid_snp_bases: ClassVar[frozenset[str]] = frozenset({"A", "C", "T", "G"})

    @computed_field
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


class CatalogScoreVariant(BaseModel):
    """A model representing a row from a PGS Catalog scoring file, defined here:
    https://www.pgscatalog.org/downloads/#scoring_columns
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
    effect_allele: Optional[Allele] = Field(
        default=None,
        title="Effect Allele",
        description="The allele that's dosage is counted (e.g. {0, 1, 2}) and multiplied by the variant's weight (effect_weight) when calculating score. The effect allele is also known as the 'risk allele'. Note: this does not necessarily need to correspond to the minor allele/alternative allele.",
    )
    other_allele: Optional[Allele] = Field(
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
    # all effect weight fields are handled as strings internally to faithfully reproduce author-uploaded scores (i.e. avoid any floating point errors)
    effect_weight: Optional[str] = Field(
        default=None,
        title="Variant Weight",
        description="Value of the effect that is multiplied by the dosage of the effect allele (effect_allele) when calculating the score. Additional information on how the effect_weight was derived is in the weight_type field of the header, and score development method in the metadata downloads.",
        coerce_numbers_to_str=True,
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
    dosage_0_weight: Optional[str] = Field(
        default=None,
        title="Effect weight with 0 copy of the effect allele",
        description="Weights that are specific to different dosages of the effect_allele (e.g. {0, 1, 2} copies) can also be reported when the the contribution of the variants to the score is not encoded as additive, dominant, or recessive. In this case three columns are added corresponding to which variant weight should be applied for each dosage, where the column name is formated as dosage_#_weight where the # sign indicates the number of effect_allele copies.",
        coerce_numbers_to_str=True,
    )
    dosage_1_weight: Optional[str] = Field(
        default=None,
        title="Effect weight with 1 copy of the effect allele",
        description="Weights that are specific to different dosages of the effect_allele (e.g. {0, 1, 2} copies) can also be reported when the the contribution of the variants to the score is not encoded as additive, dominant, or recessive. In this case three columns are added corresponding to which variant weight should be applied for each dosage, where the column name is formated as dosage_#_weight where the # sign indicates the number of effect_allele copies.",
        coerce_numbers_to_str=True,
    )
    dosage_2_weight: Optional[str] = Field(
        default=None,
        title="Effect weight with 2 copies of the effect allele",
        description="Weights that are specific to different dosages of the effect_allele (e.g. {0, 1, 2} copies) can also be reported when the the contribution of the variants to the score is not encoded as additive, dominant, or recessive. In this case three columns are added corresponding to which variant weight should be applied for each dosage, where the column name is formated as dosage_#_weight where the # sign indicates the number of effect_allele copies.",
        coerce_numbers_to_str=True,
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
    hm_inferOtherAllele: Optional[Allele] = Field(
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

    @computed_field
    @cached_property
    def variant_id(self) -> str:
        """ID = chr:pos:effect_allele:other_allele"""
        return ":".join(
            [
                str(getattr(self, k) or "")  # correctly handles None elements
                for k in ["chr_name", "chr_position", "effect_allele", "other_allele"]
            ]
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

    @field_validator(
        "effect_weight", "dosage_0_weight", "dosage_1_weight", "dosage_2_weight"
    )
    @classmethod
    def effect_weight_must_float(cls, weight: str) -> str:
        _ = float(weight)  # will raise a ValueError if conversion fails
        return weight

    @field_validator(
        "effect_allele", "other_allele", "hm_inferOtherAllele", mode="before"
    )
    @classmethod
    def alleles_must_parse(cls, value):
        if isinstance(value, str):
            return Allele(allele=value)
        else:
            raise ValueError(f"Can't parse {value=}")

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
            case str(), str(), str(), str():
                raise ValueError("Additive and non-additive fields are present")
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
                raise TypeError(
                    f"Bad position: {self.rsID=}, {self.chr_name=}, {self.chr_position=}"
                )

        return self

    @model_validator(mode="after")
    def check_harmonised_columns(self) -> Self:
        if self.is_harmonised:
            for x in self.harmonised_columns:
                if getattr(self, x) is None:
                    raise TypeError(f"Missing harmonised column data: {x}")
        return self


class ScoreFormatVersion(str, enum.Enum):
    v2 = "2.0"


class WeightType(str, enum.Enum):
    BETA = "beta"
    ODDSRATIO = "OR"
    HAZARDRATIO = "HR"
    NR = "NR"


class ScoreHeader(BaseModel):
    """Headers store useful metadata about a scoring file.

    Data validation is less strict than the CatalogScoreHeader, to make
    it easier for people to use custom scoring files with the PGS Catalog Calculator.

    >>> ScoreHeader(**{"pgs_id": "PGS123456", "trait_reported": "testtrait", "genome_build": "GRCh38"})
    ScoreHeader(pgs_id='PGS123456', pgs_name=None, trait_reported='testtrait', genome_build=GenomeBuild.GRCh38)
    """

    pgs_id: str = Field(title="PGS identifier")
    pgs_name: Optional[str] = Field(description="PGS name", default=None)
    trait_reported: str = Field(description="Trait name")
    # genome build is Optional because "NR" is represented internally as None
    genome_build: Optional[GenomeBuild] = Field(description="Genome build")

    @classmethod
    def _parse_genome_build(cls, value: str) -> Optional[GenomeBuild]:
        if value == "NR":
            return None
        else:
            return GenomeBuild.from_string(value)

    @field_validator("genome_build", mode="before")
    @classmethod
    def parse_genome_build(cls, value: str) -> Optional[GenomeBuild]:
        return cls._parse_genome_build(value)

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

            return cls(**header)


class CatalogScoreHeader(ScoreHeader):
    """A ScoreHeader that validates the PGS Catalog Scoring File header standard

    https://www.pgscatalog.org/downloads/#dl_ftp_scoring

    >>> from ._config import Config
    >>> testpath = Config.ROOT_DIR / "tests" / "data" / "PGS000001_hmPOS_GRCh38.txt.gz"
    >>> ScoreHeader.from_path(testpath) # doctest: +ELLIPSIS
    ScoreHeader(format_version=<ScoreFormatVersion.v2: '2.0'>, pgs_id='PGS000001', pgs_name='PRS77_BC', trait_reported='Breast cancer', trait_mapped=['breast carcinoma'], trait_efo=['EFO_0000305'], genome_build=None, variants_number=77, weight_type=None)
    """

    format_version: ScoreFormatVersion
    trait_mapped: list[str] = Field(description="Trait name")
    trait_efo: list[str] = Field(
        description="Ontology trait name, e.g. 'breast carcinoma"
    )
    variants_number: int = Field(
        gt=0, description="Number of variants listed in the PGS", default=None
    )
    # note: we'll make sure to serialise None values here and in genome_build as string "NR"
    weight_type: Optional[WeightType] = Field(
        description="Variant weight type", default=None
    )

    ##SOURCE INFORMATION
    pgp_id: str
    citation: str
    ##HARMONIZATION DETAILS
    HmPOS_build: Optional[GenomeBuild] = Field(default=None)
    HmPOS_date: Optional[date] = Field(default=None)
    HmPOS_match_pos: Optional[str] = Field(default=None)
    HmPOS_match_chr: Optional[str] = Field(default=None)

    # note: only included when different from default
    license: Optional[str] = Field(
        "PGS obtained from the Catalog should be cited appropriately, and "
        "used in accordance with any licensing restrictions set by the authors. See "
        "EBI Terms of Use (https://www.ebi.ac.uk/about/terms-of-use/) for additional "
        "details.",
        repr=False,
    )

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

    @field_validator("pgs_id")
    @classmethod
    def check_pgp_id(cls, pgp_id: str) -> str:
        return cls._check_accession(pgp_id, "PGP")

    @field_validator("genome_build", "HmPOS_build", mode="before")
    @classmethod
    def parse_genome_build(cls, value: str) -> Optional[GenomeBuild]:
        return cls._parse_genome_build(value)

    @field_validator("format_version")
    @classmethod
    def check_format_version(cls, version: ScoreFormatVersion) -> ScoreFormatVersion:
        if version != ScoreFormatVersion.v2:
            raise ValueError(f"Invalid format_version: {version}")
        return version

    @field_validator("weight_type")
    @classmethod
    def parse_weight_type(cls, value: WeightType) -> Optional[WeightType]:
        if value == WeightType.NR:
            value = None
        return value
