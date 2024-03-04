import logging
import pathlib
from collections import namedtuple

import pandas as pd

from .ancestry.tools import choose_pval_threshold, compare_ancestry
from .principalcomponents import PopulationType
from .ancestry import read


logger = logging.getLogger(__name__)


class AggregatedPGS:
    """A PGS that's been aggregated and melted, and may contain multiple samplesets

    >>> from ._config import Config
    >>> score_path = Config.ROOT_DIR / "tests" / "aggregated_scores.txt.gz"
    >>> AggregatedPGS(path=score_path, target_name="hgdp")
    AggregatedPGS(path=PosixPath('.../pgscatalog.calclib/tests/aggregated_scores.txt.gz'))

    """

    def __init__(self, *, target_name, df=None, path=None):
        if df is None and path is None:
            raise TypeError("df or path is required")

        self._path = path
        self._df = None
        self._target_name = target_name

    @property
    def path(self):
        return self._path

    @property
    def target_name(self):
        return self._target_name

    @property
    def df(self):
        if self._df is None:
            self._df = read.read_pgs(self._path)
        return self._df

    def __repr__(self):
        return f"{type(self).__name__}(path={repr(self.path)})"

    def _check_overlap(self, ref_pc, target_pc):
        """Before adjusting, there should be perfect target sample overlap"""
        pca_ref_samples = set(ref_pc.df.index.get_level_values(1))
        pca_target_samples = set(target_pc.df.index.get_level_values(1))
        score_ref_samples = set(self.df.loc["reference"].index)
        score_target_samples = set(self.df.loc[self.target_name].index)

        if not pca_ref_samples.issubset(score_ref_samples):
            logger.critical(
                "Error: PGS data missing for reference samples with PCA data"
            )
            raise ValueError

        if not pca_target_samples.issubset(score_target_samples):
            logger.critical("Error: PGS data missing for target samples with PCA data.")
            raise ValueError

    def adjust(self, *, ref_pc, target_pc, **kwargs):
        """
        >>> from ._config import Config
        >>> ref_pc = PrincipalComponents(pcs_path=[Config.ROOT_DIR / "tests" / "ref.pcs"], dataset="reference", psam_path=Config.ROOT_DIR / "tests" / "ref.psam", pop_type=PopulationType.REFERENCE)
        >>> target_pcs = PrincipalComponents(pcs_path=Config.ROOT_DIR / "tests" / "target.pcs", dataset="target", pop_type=PopulationType.TARGET)
        >>> score_path = Config.ROOT_DIR / "tests" / "aggregated_scores.txt.gz"
        >>> AggregatedPGS(path=score_path, target_name="hgdp").adjust(ref_pc=ref_pc, target_pc=target_pcs)
        """
        if ref_pc.pop_type != PopulationType.REFERENCE:
            raise ValueError("ref_pc argument has bad pop_type")

        if target_pc.pop_type != PopulationType.TARGET:
            raise ValueError("target_pc argument has bad pop_type")

        self._check_overlap(ref_pc=ref_pc, target_pc=target_pc)

        # join pgs + pca data
        target_df = target_pc.df.join(self.df.loc[self.target_name], on="IID")
        reference_df = ref_pc.df.join(self.df.loc["reference"], on="IID")

        # set up
        ancestry_args = namedtuple("ancestry_args", ["method_compare", "pThreshold"])
        args = ancestry_args(
            kwargs.get("method_compare", "RandomForest"), kwargs.get("pThreshold", None)
        )
        assignment_threshold_p = choose_pval_threshold(args)

        # TODO: bork
        ancestry_ref, ancestry_target, compare_info = compare_ancestry(
            ref_df=reference_df,
            ref_pop_col=ref_pc.poplabel,
            ref_train_col="Unrelated",
            target_df=target_df,
            n_pcs=ref_pc.npcs_popcomp,
            method=args.method_compare,
            p_threshold=assignment_threshold_p,
        )

        pass


class PolygenicScore:
    """Represents the output of plink2 --score written to a file

    >>> from ._config import Config
    >>> import reprlib
    >>> score1 = Config.ROOT_DIR / "tests" / "cineca_22_additive_0.sscore.zst"
    >>> pgs1 = PolygenicScore(sampleset="test", path=score1)  # doctest: +ELLIPSIS
    >>> pgs1
    PolygenicScore(sampleset='test', path=PosixPath('.../cineca_22_additive_0.sscore.zst'))
    >>> pgs2 = PolygenicScore(sampleset="test", path=score1)
    >>> reprlib.repr(pgs1.read().to_dict()) # doctest: +ELLIPSIS
    "{'DENOM': {('test', 'HG00096'): 1564, ('test', 'HG00097'): 1564, ('test', 'HG00099'): 1564, ('test', 'HG00100'): 1564, ...}, 'PGS001229_22_SUM': {('test', 'HG00096'): 0.54502, ('test', 'HG00097'): 0.674401, ('test', 'HG00099'): 0.63727, ('test', 'HG00100'): 0.863944, ...}}"

    It's often helpful to combine PGS that were split per chromosome or by effect type:

    >>> aggregated_score = pgs1 + pgs2
    >>> aggregated_score  # doctest: +ELLIPSIS
    PolygenicScore(sampleset='test', path=None)

    Once a score has been fully aggregated it can be helpful to recalculate an average:

    >>> reprlib.repr(aggregated_score.average().to_dict())  # doctest: +ELLIPSIS
    "{'DENOM': {('test', 'HG00096'): 3128, ('test', 'HG00097'): 3128, ('test', 'HG00099'): 3128, ('test', 'HG00100'): 3128, ...}, 'PGS001229_22_AVG': {('test', 'HG00096'): 0.0003484782608695652, ('test', 'HG00097'): 0.00043120268542199493, ('test', 'HG00099'): 0.0004074616368286445, ('test', 'HG00100'): 0.0005523938618925831, ...}}"

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

        for chunk in pd.read_csv(
            self.path, chunksize=self._chunksize, sep="\t", converters={"IID": str}
        ):
            df = chunk.assign(sampleset=self.sampleset).set_index(["sampleset", "#IID"])

            df.index.names = ["sampleset", "IID"]
            df = df[_select_agg_cols(df.columns)]
            yield df

    def read(self):
        """Eagerly load a PGS into a pandas dataframe"""
        if self.path is None:
            raise ValueError("Missing path")

        df = (
            pd.read_csv(self.path, sep="\t", converters={"IID": str})
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
            avgs = chunk.filter(regex="SUM$")
            avgs = avgs.divide(chunk.DENOM, axis=0)
            avgs.insert(0, "DENOM", chunk.DENOM)
            avgs.columns = avgs.columns.str.replace("_SUM", "_AVG")
            chunk_list.append(avgs)
        return pd.concat(chunk_list, axis=0)

    def melt(self):
        """Melt dataframe from wide format to long format"""
        sum_df = _melt(pd.concat(self.df, axis=0), value_name="SUM")
        avg_df = _melt(self.average(), value_name="AVG")
        return pd.concat([sum_df, avg_df.AVG], axis=1)


def _select_agg_cols(cols):
    """Select aggregatable columns"""
    keep_cols = ["DENOM"]
    return [
        x
        for x in cols
        if (x.endswith("_SUM") and (x != "NAMED_ALLELE_DOSAGE_SUM")) or (x in keep_cols)
    ]


def _melt(df, value_name):
    df = df.melt(
        id_vars=["DENOM"],
        value_name=value_name,
        var_name="PGS",
        ignore_index=False,
    )
    # e.g. PGS000822_SUM -> PGS000822
    df["PGS"] = df["PGS"].str.replace(f"_{value_name}", "")
    return df
