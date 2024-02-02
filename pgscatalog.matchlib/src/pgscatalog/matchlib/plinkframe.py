import collections.abc

from ._match.plink import plinkify


class PlinkFrame:
    """A PlinkFrame contains matches that have been formatted for plink2 --score"""

    def __init__(self, effect_type, n, df):
        self.effect_type = effect_type
        self.n = n
        self.df = df

    def __repr__(self):
        return f"{type(self).__name__}(effect_type={self.effect_type!r}, n={self.n!r}, df={self.df!r})"

    def split(self):
        return self.df.collect().partition_by("chr_name", as_dict=True)


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
