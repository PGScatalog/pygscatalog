""" This module contains classes to work with target variants. When a scoring file is
being reused to calculate scores for new genotypes, the new genotypes are target
genomes."""
import enum
import csv
import io

from pgscatalog.corelib._read import auto_open
import zstandard as zstd


class TargetVariant:
    """A single target variant, including genomic coordinates and allele information

    >>> a = TargetVariant(chrom="1", pos=12, ref="A", alt="C")
    >>> a
    TargetVariant(chrom='1', pos=12, ref='A', alt='C')
    >>> b = a
    >>> b == a
    True
    """

    def __init__(self, *, chrom, pos, ref, alt, **kwargs):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref
        self.alt = alt

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(chrom={repr(self.chrom)}, pos="
            f"{repr(self.pos)}, "
            f"ref={repr(self.ref)}, alt={repr(self.alt)})"
        )

    def __hash__(self):
        return hash((self.chrom, self.pos, self.ref, self.alt))

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

    >>> from pgscatalog.corelib import Config  # ignore, only to load test data
    >>> pvar = TargetVariants(Config.ROOT_DIR / "tests" / "hapnest.pvar")
    >>> pvar.ftype
    <TargetType.PVAR: 1>

    Variants are generators:

    >>> pvar.variants # doctest: +ELLIPSIS
    <generator object read_pvar at ...>
    >>> for variant in pvar.variants:
    ...   variant
    ...   break
    TargetVariant(chrom='14', pos=65003549, ref='T', alt='C')

    gzip and zstandard compression is transparently handled for pvar:

    >>> pvar = TargetVariants(Config.ROOT_DIR / "tests" / "hapnest.pvar.zst")
    >>> for variant in pvar.variants:
    ...   variant
    ...   break
    TargetVariant(chrom='14', pos=65003549, ref='T', alt='C')

    The same is true for bim files:

    >>> bim = TargetVariants(Config.ROOT_DIR / "tests" / "hapnest.bim.gz")
    >>> bim.ftype
    <TargetType.BIM: 2>
    >>> for variant in bim.variants:
    ...    variant
    ...    break
    TargetVariant(chrom='1', pos=10180, ref='C', alt='T')
    >>> bim = TargetVariants(Config.ROOT_DIR / "tests" / "hapnest.bim.zst")
    >>> for variant in bim.variants:
    ...    variant
    ...    break
    TargetVariant(chrom='1', pos=10180, ref='C', alt='T')
    >>> bim = TargetVariants(Config.ROOT_DIR / "tests" / "hapnest.bim")
    >>> for variant in bim.variants:
    ...    variant
    ...    break
    TargetVariant(chrom='1', pos=10180, ref='C', alt='T')

    Note, A1/A2 isn't guaranteed to be ref/alt because of PLINK1 file format
    limitations. PGS Catalog libraries handle this internally, but you should be
    aware REF/ALT can be swapped by plink during VCF to bim conversion.
    """

    def __init__(self, path):
        match n := path.name:
            case _ if "pvar" in n:
                self.ftype = TargetType.PVAR
            case _ if "bim" in n:
                self.ftype = TargetType.BIM
            case _:
                raise ValueError(f"Unknown target type {n!r}")

        self._path = str(path)

    def __repr__(self):
        return f"{self.__class__.__name__}(path={repr(self.path)})"

    @property
    def path(self):
        return self._path

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
    try:
        with auto_open(path, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            yield from yield_pvar_rows(reader)
    except UnicodeDecodeError:
        with open(path, "rb") as f:
            dctx = zstd.ZstdDecompressor()
            stream_reader = dctx.stream_reader(f)
            text_stream = io.TextIOWrapper(stream_reader, encoding="utf-8")
            reader = csv.DictReader(text_stream, delimiter="\t")
            yield from yield_pvar_rows(reader)


def read_bim(path):
    try:
        with auto_open(path, "rt") as f:
            reader = csv.reader(f, delimiter="\t")
            yield from yield_bim_rows(reader)
    except UnicodeDecodeError:
        with open(path, "rb") as f:
            dctx = zstd.ZstdDecompressor()
            stream_reader = dctx.stream_reader(f)
            text_stream = io.TextIOWrapper(stream_reader, encoding="utf-8")
            reader = csv.reader(text_stream, delimiter="\t")
            yield from yield_bim_rows(reader)


def yield_pvar_rows(iterable):
    fields = {"#CHROM": "chrom", "POS": "pos", "REF": "ref", "ALT": "alt"}
    for row in iterable:
        # update dict keys to match expected kwargs in TargetVariant
        for k, v in fields.items():
            row[v] = row.pop(k)
        yield TargetVariant(**row)


def yield_bim_rows(iterable):
    # yes, A1/A2 in bim isn't ref/alt
    fields = ["chrom", "ID", "pos_cm", "pos", "ref", "alt"]
    for row in iterable:
        yield TargetVariant(**dict(zip(fields, row, strict=True)))


class TargetType(enum.Enum):
    PVAR = enum.auto()
    BIM = enum.auto()
