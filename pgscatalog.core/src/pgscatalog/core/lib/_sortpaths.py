""" This module assumes you're working with paths that follow the format:

{sampleset}_{chrom}_{effect_type}_{n}
"""
from natsort import natsort_keygen, ns


def effect_type_keyfunc():
    """Return a key that sorts by effect type and n. Chromosome order doesn't matter.

    This is useful for things like itertools.groupby which expect sorted input

    >>> import pathlib
    >>> paths = [pathlib.Path("ukb_2_dominant_0.txt.gz"), pathlib.Path("ukb_X_additive_0.txt.gz"), pathlib.Path("ukb_X_additive_1.txt.gz"), pathlib.Path("ukb_1_recessive_0.txt.gz")]
    >>> sorted(paths, key=effect_type_keyfunc())
    [PosixPath('ukb_X_additive_0.txt.gz'), PosixPath('ukb_X_additive_1.txt.gz'), PosixPath('ukb_2_dominant_0.txt.gz'), PosixPath('ukb_1_recessive_0.txt.gz')]
    """
    return natsort_keygen(key=lambda x: x.stem.split("_")[2:], alg=ns.REAL)


def chrom_keyfunc():
    """Return a key that sorts by chromosome, including non-integer chromosomes

    >>> import pathlib
    >>> paths = [pathlib.Path("ukb_2_additive_0.txt.gz"), pathlib.Path("ukb_X_additive_0.txt.gz"), pathlib.Path("ukb_1_additive_0.txt.gz")]
    >>> sorted(paths, key=chrom_keyfunc())
    [PosixPath('ukb_1_additive_0.txt.gz'), PosixPath('ukb_2_additive_0.txt.gz'), PosixPath('ukb_X_additive_0.txt.gz')]
    """
    return natsort_keygen(key=lambda x: x.stem.split("_")[1], alg=ns.REAL)
