""" This module contains internal classes that store matched variants and format them
to be compatible with plink2 --score
"""
import collections.abc
import gzip
import pathlib

from ._match.plink import plinkify, pivot_score


class PlinkFrame:
    """A PlinkFrame contains matches and makes them suitable for plink2 --score

    A PlinkFrame represents the largest possible set of variants that can be written
    to one plink score file. Including other effect types breaks dosage calculation.
    Including duplicate IDs with different effect alleles (n) breaks plink2 --score
    When writing to a file the set of variants can be split into separate files (one
    per chromosome) for larger datasets.
    """

    def __init__(self, effect_type, n, df):
        self.effect_type = effect_type
        self.n = n
        self.df = df

    def __repr__(self):
        return f"{type(self).__name__}(effect_type={self.effect_type!r}, n={self.n!r}, df={self.df!r})"

    def split_pivot(self):
        """Splitting scoring files is helpful for split - apply - combine on big data"""
        dfs = self.df.collect().partition_by(["chr_name"], as_dict=True)
        return {k[0]: v.pipe(pivot_score) for k, v in dfs.items()}

    def pivot_wide(self):
        """Pivoting wide is important to enable parallel score calculation"""
        return self.df.collect().pipe(pivot_score)

    def write(self, directory, dataset, split=False):
        """Write a plink2 --score compatible file, optionally splitting"""
        fouts = []
        if split:
            dfs = self.split_pivot()
            for chrom, df in dfs.items():
                fout = (
                    pathlib.Path(directory)
                    / f"{dataset}_{chrom}_{str(self.effect_type)}_{self.n}.scorefile.gz"
                )
                with gzip.open(fout, "wb") as f:
                    df.write_csv(f, separator="\t")
                fouts.append(fout)
        else:
            chrom = "ALL"
            fout = (
                pathlib.Path(directory)
                / f"{dataset}_{chrom}_{str(self.effect_type)}_{self.n}.scorefile.gz"
            )
            df = self.pivot_wide()
            with gzip.open(fout, "wb") as f:
                df.write_csv(f, separator="\t")
            fouts.append(fout)
        return fouts


class PlinkFrames(collections.abc.Sequence):
    """A container of PlinkFrames"""

    def __init__(self, *elements):
        self._elements = list(elements)

    def __getitem__(self, item):
        return self._elements[item]

    def __len__(self):
        return len(self._elements)

    def __repr__(self):
        return f"{type(self).__name__}({self._elements!r})"

    @classmethod
    def from_matchresult(cls, matchresult):
        effect_types = plinkify(matchresult)
        plinkframes = []
        for effect_type, dataframes in effect_types.items():
            for i, df in enumerate(dataframes):
                plinkframes.append(PlinkFrame(effect_type=effect_type, n=i, df=df))

        return cls(*plinkframes)
