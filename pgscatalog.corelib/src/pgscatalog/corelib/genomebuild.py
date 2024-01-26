import enum


class GenomeBuild(enum.Enum):
    """Enumeration of genome build: the reference genome release that a scoring file
    is aligned to.

    >>> GenomeBuild.GRCh38
    GenomeBuild.GRCh38
    """

    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"
    # just included to handle older files, incompatible unless harmonised:
    NCBI36 = "NCBI36"  # ew

    def __str__(self):
        return str(self.value)

    def __repr__(self):
        return f"{type(self).__name__}.{self.name}"

    @classmethod
    def from_string(cls, build):
        """
        :param build: genome build string
        :return: :class:`GenomeBuild`
        :raises ValueError: From an unsupported build string

        >>> GenomeBuild.from_string("GRCh38")
        GenomeBuild.GRCh38
        >>> str(GenomeBuild.from_string("GRCh37"))
        'GRCh37'
        >>> GenomeBuild.from_string("NR") is None
        True
        >>> GenomeBuild.from_string("pangenome")
        Traceback (most recent call last):
        ...
        ValueError: Can't match build='pangenome'
        """
        match build:
            case "GRCh37" | "hg19":
                return cls(GenomeBuild.GRCh37)
            case "GRCh38" | "hg38":
                return cls(GenomeBuild.GRCh38)
            case "NR" | "" | None:
                return None
            case "NCBI36" | "hg18":
                return cls(GenomeBuild.NCBI36)
            case _:
                raise ValueError(f"Can't match {build=}")
