"""This module contains classes to work with target variants. When a scoring file is
being reused to calculate scores for new genotypes, the new genotypes are target
genomes."""

import enum
import csv
import pathlib

from xopen import xopen


class TargetVariant:
    """A single target variant, including genomic coordinates and allele information

    >>> a = TargetVariant(chrom="1", pos=12, ref="A", alt="C", id='1:12:A:C')
    >>> a
    TargetVariant(chrom='1', pos=12, ref='A', alt='C', id='1:12:A:C')
    >>> b = a
    >>> b == a
    True
    """

    def __init__(self, *, chrom, pos, ref, alt, id):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref
        self.alt = alt
        self.id = id

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(chrom={repr(self.chrom)}, pos="
            f"{repr(self.pos)}, "
            f"ref={repr(self.ref)}, "
            f"alt={repr(self.alt)}, "
            f"id={repr(self.id)})"
        )

    def __hash__(self):
        return hash((self.chrom, self.pos, self.ref, self.alt, self.id))

    def __eq__(self, other):
        if isinstance(other, TargetVariant):
            return (self.chrom, self.pos, self.ref, self.alt) == (
                other.chrom,
                other.pos,
                other.ref,
                other.alt,
            )
        return False


class TargetVariants:
    """A container of :class:`TargetVariant`

    :raises: FileNotFoundError

    >>> from pgscatalog.core import Config  # ignore, only to load test data
    >>> pvar = TargetVariants(Config.ROOT_DIR / "tests" / "data" / "hapnest.pvar")
    >>> pvar.ftype
    <TargetType.PVAR: 1>

    Iterating over TargetVariants is done via a read-only generator attribute:

    >>> pvar.variants # doctest: +ELLIPSIS
    <generator object read_pvar at ...>
    >>> for variant in pvar:
    ...   variant
    ...   break
    TargetVariant(chrom='14', pos=65003549, ref='T', alt='C', id='14:65003549:T:C')

    gzip and zstandard compression is transparently handled for pvar:

    >>> pvar = TargetVariants(Config.ROOT_DIR / "tests" / "data" / "hapnest.pvar.zst")
    >>> for variant in pvar:
    ...   variant
    ...   break
    TargetVariant(chrom='14', pos=65003549, ref='T', alt='C', id='14:65003549:T:C')

    The same is true for bim files:

    >>> bim = TargetVariants(Config.ROOT_DIR / "tests" / "data" / "hapnest.bim.gz")
    >>> bim.ftype
    <TargetType.BIM: 2>
    >>> for variant in bim:
    ...    variant
    ...    break
    TargetVariant(chrom='1', pos=10180, ref='C', alt='T', id='1:10180:T:C')
    >>> bim = TargetVariants(Config.ROOT_DIR / "tests" / "data" / "hapnest.bim.zst")
    >>> for variant in bim:
    ...    variant
    ...    break
    TargetVariant(chrom='1', pos=10180, ref='C', alt='T', id='1:10180:T:C')
    >>> bim = TargetVariants(Config.ROOT_DIR / "tests" / "data" / "hapnest.bim")
    >>> for variant in bim:
    ...    variant
    ...    break
    TargetVariant(chrom='1', pos=10180, ref='C', alt='T', id='1:10180:T:C')

    Note, A1/A2 isn't guaranteed to be ref/alt because of PLINK1 file format
    limitations. PGS Catalog libraries handle this internally, but you should be
    aware REF/ALT can be swapped by plink during VCF to bim conversion.

    Some pvar files can contain a lot of comments in the header, which are ignored:

    >>> pvar = TargetVariants(Config.ROOT_DIR / "tests" / "data" / "1000G.pvar")
    >>> for variant in pvar:
    ...     variant
    ...     break
    TargetVariant(chrom='1', pos=10390, ref='CCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA', alt='C', id='1:10390:CCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA:C')
    """

    def __init__(self, path, chrom=None):
        match n := pathlib.Path(path).name:
            case _ if "pvar" in n:
                self.ftype = TargetType.PVAR
            case _ if "bim" in n:
                self.ftype = TargetType.BIM
            case _:
                raise ValueError(f"Unknown target type {n!r}")

        self._chrom = chrom
        self._path = str(path)

    def __repr__(self):
        return f"{self.__class__.__name__}(path={repr(self.path)})"

    def __iter__(self):
        yield from self.variants

    @property
    def path(self):
        return self._path

    @property
    def chrom(self):
        return self._chrom

    @property
    def variants(self):
        match self.ftype:
            case TargetType.BIM:
                return read_bim(self.path)
            case TargetType.PVAR:
                return read_pvar(self.path)
            case _:
                raise ValueError


def read_pvar(path):
    """Read plink2 pvar variant information files using python core library"""

    with xopen(path, "rt") as f:
        # pvars do have a header column and support arbitrary columns
        for line in f:
            if line.startswith("##"):
                continue
            else:
                fieldnames = line.strip().split("\t")
                break
        reader = csv.DictReader(f, fieldnames=fieldnames, delimiter="\t")
        fields = {
            "#CHROM": "chrom",
            "POS": "pos",
            "REF": "ref",
            "ALT": "alt",
            "ID": "id",
        }
        for row in reader:
            yield TargetVariant(**{v: row[k] for k, v in fields.items()})


def read_bim(path):
    """Read plink1 bim variant information files using python core library"""
    with xopen(path, "rt") as f:
        # bims don't have header column
        reader = csv.reader(f, delimiter="\t")
        # yes, A1/A2 in bim isn't ref/alt
        fields = ["chrom", "id", "pos_cm", "pos", "ref", "alt"]
        for row in reader:
            row = dict(zip(fields, row, strict=True))
            yield TargetVariant(
                **{k: row[k] for k in ("chrom", "pos", "ref", "alt", "id")}
            )


class TargetType(enum.Enum):
    PVAR = enum.auto()
    BIM = enum.auto()
