from enum import Enum


class EffectType(Enum):
    """Enums that represent inheritance models. The vast majority of variants have
    an additive effect type.

    This changes downstream PGS calculation:

    * ScoreVariants with an additive effect type will always be added to the PGS sum.
    * ScoreVariants with a dominant effect type are only added to the PGS sum if there is at least one copy of the effect allele.
    * ScoreVariants with a recessive effect type are only added to the PGS sum if there are two copies of the effect allele.

    Non-additive variants aren't supported and will raise an exception

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
