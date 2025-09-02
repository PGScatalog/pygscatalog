import re
from typing import ClassVar, Any, Optional

from pydantic import (
    BaseModel,
    model_validator,
    ValidationError, field_validator
)
from typing_extensions import Self

from pgscatalog.core.lib.models import ScoreVariant
import pgscatalog.validate.lib.errors as errors


class ScoringFileValidationError(Exception):
    """Wrapper class for pydantic validation errors that also includes the row number where the error occurred."""

    row: int
    error: ValidationError

    def __init__(self, row, error: ValidationError):
        self.row = row
        self.error = error


class ColumnNames(BaseModel):
    """
    Class for validating the column names of a scoring file. This is an important step prior to
    validating all the variants as a simply misnamed column needs to be reported only once
    rather than for each variant.
    """

    columns: set[str]
    hm: bool = False

    base_columns: ClassVar[set[str]] = {
        'rsID',
        'chr_name',
        'chr_position',
        'effect_allele',  # required all times
        'effect_weight',
        'other_allele',
        'variant_description',
        'weight_type',  # allowed in raw files only
        'variant_type',
        'imputation_method',
        'inclusion_criteria',
        'is_diplotype',
        'is_haplotype',
        'locus_name',
        'dosage_0_weight',
        'dosage_1_weight',
        'dosage_2_weight',
        'is_interaction',
        'is_dominant',
        'is_recessive',
        'OR',
        'HR',
        'allelefrequency_effect'
    }

    hm_columns: ClassVar[set[str]] = {
        'hm_source',
        'hm_rsID',
        'hm_chr',
        'hm_pos',
        'hm_inferOtherAllele'
    }

    @model_validator(mode='after')
    def validate_columns(self) -> Self:
        columns = self.columns

        # rsID or (chr_pos and chr_name)
        if not ('rsID' in columns
                or {'chr_name', 'chr_position'}.issubset(columns)
                or {'chr_name', 'variant_type'}.issubset(columns)  # If only complex alleles
                ):
            raise ValueError(errors.MISSING_COLUMN_RSID_OR_COORD)

        # effect_weight or dosage_{0..2}_weight
        if not ('effect_weight' in columns or
                {'dosage_0_weight', 'dosage_1_weight', 'dosage_2_weight'}.issubset(columns)):
            raise ValueError(errors.MISSING_COLUMN_EFFECT_OR_DOSAGE_WEIGHT)

        # allelefrequency_effect_[Ancestry] variable columns
        pattern = r'^allelefrequency_effect_'
        columns = {c for c in columns if not re.match(pattern, c)}

        # Don't allow non-listed columns
        allowed_columns = self.base_columns.union(self.hm_columns) if self.hm else self.base_columns
        if not columns.issubset(allowed_columns):
            raise ValueError(errors.UNEXPECTED_COLUMNS.format(columns=str(columns - set(allowed_columns))))

        return self


class ValidationVariant(ScoreVariant):
    """Adds an extra layer of validation of scoring file variants"""

    ## Allowed chromosome names
    # Either 1..22, X, Y, XY or MT
    allowed_chr_names: ClassVar[set[str]] = {
        *map(str, range(1, 23)),
        'X',
        'Y',
        'XY',
        'MT'
    }
    # Or random/alt contigs like "1_KI270766v1_alt"
    _contig_pattern: ClassVar[re.Pattern] = re.compile(
        rf"^({'|'.join(map(re.escape, allowed_chr_names))})_([^_]+)_(alt|random)$"
    )

    @staticmethod
    def __value_contains_non_printable_characters(value) -> bool:
        return any([not c.isprintable() for c in str(value)])

    @staticmethod
    def __value_contains_leading_or_trailing_spaces(value) -> bool:
        return str(value).strip() != str(value)

    @field_validator(
        *ColumnNames.base_columns, *ColumnNames.hm_columns,
        mode="before",
        check_fields=False,
    )
    @classmethod
    def extra_field_validation(cls, v: Any) -> Optional[Any]:
        """
        Validate all fields for common mistakes such as non-printable characters or leading/trailing spaces.
        """
        if cls.__value_contains_non_printable_characters(v):
            raise ValueError(errors.NON_PRINTABLE_CHAR.format(input=repr(v)))
        if cls.__value_contains_leading_or_trailing_spaces(v):
            raise ValueError(errors.LEADING_OR_TRAILING_SPACE.format(input=repr(v)))
        return v

    @field_validator('chr_name', 'hm_chr', mode="after")
    @classmethod
    def validate_chr_name(cls, chr_name: str) -> str:
        """
        Ensuring that chromosome names are valid (1..22, X, Y, XY, MT or a valid contig name)
        """
        if chr_name.upper() not in cls.allowed_chr_names:
            # Allow alternate/random contigs like "1_KI270766v1_alt". Might be found in harmonised files
            if not cls._contig_pattern.match(chr_name):
                raise ValueError(errors.INVALID_CHR_NAME.format(input=chr_name))
        return chr_name


class ComplexVariant(ValidationVariant):
    """Class for validating complex variants such as HLA, CYP and APOE variants."""

    @model_validator(mode='after')
    def validate(self) -> Self:

        match self.variant_type:
            case "APOE_allele":
                self._validate_apoe()
            case v if v.startswith('HLA_'):
                self._validate_hla()
            case "CYP_allele":
                self._validate_cyp()
            case _:
                raise ValueError(errors.INVALID_VARIANT_TYPE.format(input=self.variant_type))

        # Diplotype/haplotype
        if self.is_diplotype and self.is_haplotype:
            raise ValueError(errors.BOTH_DIPLOTYPE_AND_HAPLOTYPE)
        if self.is_diplotype and self.effect_allele:
            if len(self.effect_allele.allele.split("/")) != 2:
                raise ValueError(errors.MISSING_DIPLOTYPE_ALLELE)

        return self

    def _validate_hla(self):
        variant_type = self.variant_type
        if variant_type.endswith("_allele") or variant_type.endswith("_serotype"):
            # Validation of potential exclusion groups
            var_desc_exclusion_groups = self._get_exclusion_groups()
            alleles_exclusion_groups = set()
            effect_allele = self.effect_allele.allele
            for allele in re.split(r"[/;]", effect_allele):
                for field_value in re.split(r'[*:]', allele):
                    if field_value.startswith('X'):
                        alleles_exclusion_groups.add(field_value)
                        if field_value not in var_desc_exclusion_groups:
                            raise ValueError(errors.EXCLUSION_GROUP_NOT_IN_VAR_DESC.format(input=field_value))
            if set(var_desc_exclusion_groups.keys()) != alleles_exclusion_groups:
                raise ValueError(errors.EXCLUSION_GROUPS_NOT_IN_EFFECT_ALLELE
                                 .format(input=set(var_desc_exclusion_groups.keys()) - alleles_exclusion_groups))

    def _validate_apoe(self):
        if self.locus_name != 'APOE':
            raise ValueError(errors.LOCUS_NAME_NOT_APOE)
        pass

    def _validate_cyp(self):
        pass

    def _get_exclusion_groups(self) -> dict:
        exclusion_groups = {}
        variant_description_items = [p.strip() for p in re.split(r';(?![^\[]*])', self.variant_description)]
        for item in variant_description_items:
            key, value = item.split("=", 1)
            if re.match(r'X[0-9]*!', key):
                exclusion_group_id = key[:-1]
                if exclusion_group_id in exclusion_groups:
                    raise ValueError(errors.DUPLICATE_EXCLUSION_GROUPS.format(input=exclusion_group_id))
                exclusion_groups[exclusion_group_id] = value.strip('[').strip(']')
        return exclusion_groups
