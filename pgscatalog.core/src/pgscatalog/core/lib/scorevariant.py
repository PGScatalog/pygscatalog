"""This module contains classes that compose a ScoreVariant: a variant in a PGS
Catalog Scoring File."""
from enum import Enum
from typing import Optional


class EffectAllele:
    """A class that represents an effect allele found in PGS Catalog scoring files

    The allele that's dosage is counted (e.g. {0, 1, 2}) and multiplied by the variant's
    weight (effect_weight) when calculating score. The effect allele is also known as
    the 'risk allele'.

    >>> simple_ea = EffectAllele("A")
    >>> simple_ea
    EffectAllele("A")
    >>> simple_ea.is_snp
    True
    >>> str(simple_ea)
    'A'
    >>> EffectAllele("AG")
    EffectAllele("AG")
    >>> hla_example = EffectAllele("+")
    >>> hla_example
    EffectAllele("+")
    >>> hla_example.is_snp
    False
    """

    _valid_snp_bases = frozenset({"A", "C", "T", "G"})
    __slots__ = ("_allele", "_is_snp")

    def __init__(self, allele):
        self._allele = str(allele)
        self._is_snp = None  # computed when accessed

    def __repr__(self):
        return f'{type(self).__name__}("{self.allele}")'

    def __str__(self):
        return self.allele

    @property
    def allele(self):
        return self._allele

    @allele.setter
    def allele(self, value):
        self._allele = str(value)
        self._is_snp = None  # reset _is_snp when allele is changed

    @property
    def is_snp(self) -> bool:
        """SNPs are the most common type of effect allele in PGS Catalog scoring
        files. More complex effect alleles, like HLAs or APOE genes, often require
        extra work to represent in genomes. Users should be warned about complex
        effect alleles.

        >>> ea = EffectAllele("+")
        >>> ea.is_snp
        False
        >>> ea.allele = "A"
        >>> ea.is_snp
        True
        """
        if self._is_snp is None:
            self._is_snp = not frozenset(self.allele) - self._valid_snp_bases
        return self._is_snp

    def __eq__(self, other):
        if isinstance(other, EffectAllele):
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

    def __str__(self):
        return str(self.value)

    def __repr__(self):
        # pasting __repr__ output should be sufficient to construct the class
        return f"{type(self).__name__}.{self.name}"


class ScoreVariant:
    """Represents a single row in a PGS Catalog scoring file.

    It's rare to instantiate this class directly. Instead, create a
    class:`ScoringFile`  from a path and you can lazily iterate over variants.
    """

    mandatory_fields: tuple[str] = (
        "effect_allele",
        "effect_weight",
        "accession",
        "row_nr",
    )
    optional_fields: tuple[str] = (
        "chr_name",
        "chr_position",
        "rsID",
        "other_allele",
        "hm_chr",
        "hm_pos",
        "hm_inferOtherAllele",
        "hm_source",
        "is_dominant",
        "is_recessive",
        "hm_rsID",
        "hm_match_chr",
        "hm_match_pos",
        "is_duplicated",
        "effect_type",
    )
    complex_fields: tuple[str] = ("is_haplotype", "is_diplotype", "is_interaction")

    # column names for output are used by __iter__ and when writing out
    output_fields: tuple[str] = (
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

    # slots uses magic to improve speed and memory when making millions of objects
    __slots__ = mandatory_fields + optional_fields + ("is_complex",)

    def __init__(
        self,
        *,
        effect_allele: str,
        effect_weight: str,
        accession: str,
        row_nr: int,
        chr_name: str = None,
        chr_position: int = None,
        rsID: str = None,
        other_allele: str = None,
        hm_chr: str = None,
        hm_pos: int = None,
        hm_inferOtherAllele: str = None,
        hm_source: str = None,
        is_dominant: str = None,
        is_recessive: str = None,
        hm_rsID: str = None,
        hm_match_chr: str = None,
        hm_match_pos: str = None,
        is_duplicated: bool = False,
        effect_type: EffectType = EffectType.ADDITIVE,
        is_complex: bool = False,
        **kwargs,
    ):
        # start with mandatory attributes
        self.effect_allele: EffectAllele = EffectAllele(effect_allele)
        self.effect_weight: str = effect_weight
        self.accession: str = accession
        self.row_nr: int = int(row_nr)

        # now set optional fields
        self.chr_name: Optional[str] = chr_name

        # casting to int is important for arrow export
        try:
            self.chr_position: Optional[int] = int(chr_position)
        except (ValueError, TypeError):
            self.chr_position = None

        self.rsID: Optional[str] = rsID
        self.other_allele: Optional[str] = other_allele
        self.hm_chr: Optional[str] = hm_chr

        # casting to int is important when harmonised data may replace chr_position
        try:
            self.hm_pos: Optional[int] = int(hm_pos)
        except (ValueError, TypeError):
            self.hm_pos = None

        self.hm_inferOtherAllele: Optional[str] = hm_inferOtherAllele
        self.hm_source: Optional[str] = hm_source

        match is_dominant:
            case True | "True":
                self.is_dominant = True
            case None:
                self.is_dominant = None
            case _:
                self.is_dominant = False

        match is_recessive:
            case True | "True":
                self.is_recessive = True
            case None:
                self.is_recessive = None
            case _:
                self.is_recessive = False

        self.hm_rsID: Optional[str] = hm_rsID
        self.hm_match_chr: Optional[str] = hm_match_chr
        self.hm_match_pos: Optional[str] = hm_match_pos
        self.is_duplicated: Optional[bool] = is_duplicated
        self.effect_type: EffectType = effect_type

        # these fields are important to check if variants are complex
        if any([x in kwargs for x in self.complex_fields]):
            is_complex = True
        self.is_complex: bool = is_complex

    def __repr__(self):
        class_name = type(self).__name__
        values = {}

        for key in ScoreVariant.__slots__:
            values[key] = getattr(self, key, None)

        # extract str parameter for effect allele
        values["effect_allele"] = values["effect_allele"].allele

        params = ",".join([f"{k}={repr(v)}" for k, v in values.items()])
        return f"{class_name}({params})"

    def __iter__(self):
        for attr in self.output_fields:
            yield getattr(self, attr)
