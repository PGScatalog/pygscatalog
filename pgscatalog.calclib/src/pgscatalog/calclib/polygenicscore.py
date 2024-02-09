import pathlib

import pandas as pd
import reprlib


class PolygenicScore:
    """Represents the output of plink2 --score written to a file

    >>> from ._config import Config
    >>> score1 = Config.ROOT_DIR / "tests" / "cineca_22_additive_0.sscore.zst"
    >>> pgs1 = PolygenicScore(sampleset="test", path=score1)  # doctest: +ELLIPSIS
    >>> pgs1
    PolygenicScore(sampleset='test', path=PosixPath('.../cineca_22_additive_0.sscore.zst'), df=None)
    >>> pgs2 = PolygenicScore(sampleset="test", path=score1)
    >>> pgs1.read().to_dict()  # doctest: +ELLIPSIS
    {'DENOM': ...}, 'PGS001229_22_SUM': {('test', 'HG00096'): 0.54502, ('test', 'HG00097'): 0.674401, ('test', 'HG00099'): 0.63727, ('test', 'HG00100'): 0.863944, ...}}

    It's often helpful to combine PGS that were split per chromosome or by effect type:

    >>> aggregated_score = pgs1 + pgs2
    >>> aggregated_score  # doctest: +ELLIPSIS
    PolygenicScore(sampleset='test', path=None, df={'DENOM': ...}, 'PGS001229_22_SUM': {('test', 'HG00096'): 1.09004, ('test', 'HG00097'): 1.348802, ('test', 'HG00099'): 1.27454, ('test', 'HG00100'): 1.727888, ...}})

    Once a score has been fully aggregated it can be helpful to recalculate an average:

    >>> aggregated_score.average().to_dict()  # doctest: +ELLIPSIS
    {'DENOM': ...}, 'PGS001229_22_SUM': {('test', 'HG00096'): 1.09004, ...}, 'PGS001229_22_AVG': {('test', 'HG00096'): 0.000348...

    Scores can be written to a TSV file:

    >>> import tempfile, os
    >>> outd = tempfile.mkdtemp()
    >>> aggregated_score.write(str(outd))
    >>> os.listdir(outd)
    ['aggregated_scores.txt.gz']

    With support for splitting output files by sampleset:

    >>> splitoutd = tempfile.mkdtemp()
    >>> aggregated_score.write(splitoutd, split=True)
    >>> sorted(os.listdir(splitoutd), key = lambda x: x.split("_")[0])
    ['test_pgs.txt.gz']
    """

    def __init__(self, *, sampleset, path=None, df=None):
        match (path, df):
            case (None, None):
                raise ValueError("init with path or df")
            case _ if path is not None and df is not None:
                raise ValueError("choose one to init: path or df")
            case _:
                pass

        self.path = path
        self.df = df
        self.sampleset = sampleset

    def __repr__(self):
        if self.df is not None:
            df = reprlib.repr(self.df.to_dict())
        else:
            df = reprlib.repr(None)

        return f"{type(self).__name__}(sampleset={repr(self.sampleset)}, path={repr(self.path)}, df={df})"

    def read(self):
        if self.df is None:
            df = (
                pd.read_table(self.path)
                .assign(sampleset=self.sampleset)
                .set_index(["sampleset", "#IID"])
            )

            df.index.names = ["sampleset", "IID"]
            df = df[_select_agg_cols(df.columns)]
            self.df = df
        return self.df

    def average(self):
        avgs = self.df.loc[:, self.df.columns.str.endswith("_SUM")].divide(
            self.df["DENOM"], axis=0
        )
        avgs.columns = avgs.columns.str.replace("_SUM", "_AVG")
        self.df = pd.concat([self.df, avgs], axis=1)
        return self.df

    def write(self, outdir, split=False):
        outdir = pathlib.Path(outdir)
        if split:
            for sampleset, group in self.df.groupby("sampleset"):
                fout = outdir / f"{sampleset}_pgs.txt.gz"
                group.to_csv(fout, sep="\t", compression="gzip")
        else:
            fout = outdir / "aggregated_scores.txt.gz"
            self.df.to_csv(fout, sep="\t", compression="gzip")

    def __add__(self, other):
        if isinstance(other, PolygenicScore):
            if self.sampleset != other.sampleset:
                raise ValueError("Can't add PolygenicScore with different samplesets")

            df1 = self.read()
            df2 = other.read()
            sumdf = df1.add(df2, fill_value=0)
            return PolygenicScore(sampleset=self.sampleset, df=sumdf)
        else:
            return NotImplemented


def _select_agg_cols(cols):
    """Select aggregatable columns"""
    keep_cols = ["DENOM"]
    return [
        x
        for x in cols
        if (x.endswith("_SUM") and (x != "NAMED_ALLELE_DOSAGE_SUM")) or (x in keep_cols)
    ]
