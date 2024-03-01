import pathlib

import pandas as pd


class PolygenicScore:
    """Represents the output of plink2 --score written to a file

    >>> from ._config import Config
    >>> score1 = Config.ROOT_DIR / "tests" / "cineca_22_additive_0.sscore.zst"
    >>> pgs1 = PolygenicScore(sampleset="test", path=score1)  # doctest: +ELLIPSIS
    >>> pgs1
    PolygenicScore(sampleset='test', path=PosixPath('.../cineca_22_additive_0.sscore.zst'))
    >>> pgs2 = PolygenicScore(sampleset="test", path=score1)
    >>> pgs1.read().to_dict()  # doctest: +ELLIPSIS
    {'DENOM': ...}, 'PGS001229_22_SUM': {('test', 'HG00096'): 0.54502, ('test', 'HG00097'): 0.674401, ('test', 'HG00099'): 0.63727, ('test', 'HG00100'): 0.863944, ...}}

    It's often helpful to combine PGS that were split per chromosome or by effect type:

    >>> aggregated_score = pgs1 + pgs2
    >>> aggregated_score  # doctest: +ELLIPSIS
    PolygenicScore(sampleset='test', path=None)

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

    def __init__(self, *, path=None, df=None, sampleset=None):
        match (path, df):
            case (None, None):
                raise ValueError("init with path or df")
            case _ if path is not None and df is not None:
                raise ValueError("choose one to init: path or df")
            case _:
                pass

        try:
            self.path = pathlib.Path(path)
        except TypeError:
            self.path = None

        if sampleset is None:
            self.sampleset = path.stem.split("_")[0]
        else:
            self.sampleset = sampleset

        self._chunksize = 50000

        if df is not None:
            # big df is an in-memory pandas df
            self._bigdf = df
        else:
            self._bigdf = None
            self._df = None

    def __repr__(self):
        return f"{type(self).__name__}(sampleset={repr(self.sampleset)}, path={repr(self.path)})"

    def __iter__(self):
        yield from self.df

    def __add__(self, other):
        if isinstance(other, PolygenicScore):
            dfs = []
            for df1, df2 in zip(self, other, strict=True):
                sumdf = df1.add(df2, fill_value=0)
                dfs.append(sumdf)
            return PolygenicScore(sampleset=self.sampleset, df=pd.concat(dfs, axis=0))
        else:
            return NotImplemented

    @property
    def df(self):
        if self.path is not None:
            self._df = self.lazy_read()
        elif self._bigdf is not None:
            # a fake generator
            self._df = (x for x in [self._bigdf])
        return self._df

    def lazy_read(self):
        """Lazily read a PGS in chunks"""
        if self.path is None:
            raise ValueError("Missing path")

        for chunk in pd.read_csv(self.path, chunksize=self._chunksize, sep="\t"):
            df = chunk.assign(sampleset=self.sampleset).set_index(["sampleset", "#IID"])

            df.index.names = ["sampleset", "IID"]
            df = df[_select_agg_cols(df.columns)]
            yield df

    def read(self):
        """Eagerly load a PGS into a pandas dataframe"""
        if self.path is None:
            raise ValueError("Missing path")

        df = (
            pd.read_csv(self.path, sep="\t")
            .assign(sampleset=self.sampleset)
            .set_index(["sampleset", "#IID"])
        )

        df.index.names = ["sampleset", "IID"]
        df = df[_select_agg_cols(df.columns)]
        return df

    def write(self, outdir, split=False, melt=True):
        """Write PGS to a compressed TSV"""
        outdir = pathlib.Path(outdir)
        for chunk in self:
            if melt:
                chunk = _melt(chunk, "SUM")

            if split:
                for sampleset, group in chunk.groupby("sampleset"):
                    fout = outdir / f"{sampleset}_pgs.txt.gz"
                    chunk.to_csv(fout, sep="\t", compression="gzip", mode="a")
            else:
                fout = outdir / "aggregated_scores.txt.gz"
                chunk.to_csv(fout, sep="\t", compression="gzip", mode="a")

    def average(self):
        """Recalculate average.

        This is an eager operation, and immediately returns a dataframe
        """
        chunk_list = []
        for chunk in self:
            avgs = chunk.loc[:, chunk.columns.str.endswith("_SUM")].divide(
                chunk["DENOM"], axis=0
            )
            avgs.columns = avgs.columns.str.replace("_SUM", "_AVG")
            chunk_list.append(pd.concat([chunk, avgs], axis=1))
        return pd.concat(chunk_list, axis=0)


def _select_agg_cols(cols):
    """Select aggregatable columns"""
    keep_cols = ["DENOM"]
    return [
        x
        for x in cols
        if (x.endswith("_SUM") and (x != "NAMED_ALLELE_DOSAGE_SUM")) or (x in keep_cols)
    ]


def _melt(df, value_name):
    """Melt the score dataframe from wide format to long format"""
    df = df.melt(
        id_vars=["DENOM"],
        value_name=value_name,
        var_name="PGS",
        ignore_index=False,
    )
    # e.g. PGS000822_SUM -> PGS000822
    df["PGS"] = df["PGS"].str.replace(f"_{value_name}", "")
    return df
