import dataclasses
import gzip
import logging
import pathlib
import typing

import pandas as pd

from ._ancestry.tools import (
    compare_ancestry,
    choose_pval_threshold,
    pgs_adjust,
    write_model,
)
from .principalcomponents import PopulationType
from ._ancestry import read


logger = logging.getLogger(__name__)


@dataclasses.dataclass(frozen=True)
class AdjustArguments:
    """Arguments that control genetic similarity estimation and PGS adjustment

    >>> AdjustArguments(method_compare="Mahalanobis", pThreshold=None, method_normalization=("empirical", "mean"))
    AdjustArguments(method_compare='Mahalanobis', pThreshold=None, method_normalization=('empirical', 'mean'))
    """

    method_compare: str = "RandomForest"
    pThreshold: typing.Optional[float] = None
    method_normalization: tuple[str, ...] = ("empirical", "mean", "mean+var")

    def __post_init__(self):
        if self.method_compare not in {"Mahalanobis", "RandomForest"}:
            raise ValueError(f"Bad method: {self.method_compare}")

        if not set(self.method_normalization) <= {"empirical", "mean", "mean+var"}:
            raise ValueError(f"Bad normalisation: {self.method_normalization}")


@dataclasses.dataclass(frozen=True)
class AdjustResults:
    """Results returned by :class:`AggregatedPGS.adjust()`"""

    target_label: str
    models: pd.DataFrame
    model_meta: dict
    pca: pd.DataFrame
    pgs: pd.DataFrame
    scorecols: list[str]

    def write(self, directory):
        """Write model, PGS, and PCA data to a directory"""
        self._write_model(directory)
        self._write_pgs(directory)
        self._write_pca(directory)

    def _write_pca(self, directory):
        directory = pathlib.Path(directory)
        loc_popsim_out = directory / f"{self.target_label}_popsimilarity.txt.gz"
        logger.debug(f"Writing PCA and popsim results to: {loc_popsim_out}")
        self.pca.drop(self.scorecols, axis=1).to_csv(loc_popsim_out, sep="\t")

    def _write_model(self, directory):
        """Write results to a directory"""
        directory = pathlib.Path(directory)
        write_model(
            {"pgs": self.models, "compare_pcs": self.model_meta},
            str(directory / f"{self.target_label}_info.json.gz"),
        )

    def _write_pgs(self, directory):
        scorecols = self.scorecols
        directory = pathlib.Path(directory)
        loc_pgs_out = str(directory / f"{self.target_label}_pgs.txt.gz")
        with gzip.open(loc_pgs_out, "wt") as outf:
            logger.debug(
                "Writing adjusted PGS values (long format) to: {}".format(loc_pgs_out)
            )
            for i, pgs_id in enumerate(scorecols):
                df_pgs = self.pgs.loc[:, self.pgs.columns.str.endswith(pgs_id)].melt(
                    ignore_index=False
                )  # filter to 1 PGS
                df_pgs[["method", "PGS"]] = df_pgs.variable.str.split("|", expand=True)
                df_pgs = (
                    df_pgs.drop("variable", axis=1)
                    .reset_index()
                    .pivot(
                        index=["sampleset", "IID", "PGS"],
                        columns="method",
                        values="value",
                    )
                )
                if i == 0:
                    logger.debug(
                        "{}/{}: Writing {}".format(i + 1, len(scorecols), pgs_id)
                    )
                    colorder = list(df_pgs.columns)  # to ensure sort order
                    df_pgs.to_csv(outf, sep="\t")
                else:
                    logger.debug(
                        "{}/{}: Appending {}".format(i + 1, len(scorecols), pgs_id)
                    )
                    df_pgs[colorder].to_csv(outf, sep="\t", header=False)


class AggregatedPGS:
    """A PGS that's been aggregated, melted, and probably contains samples from a reference panel and a target population.

    The most useful method in this class adjusts PGS based on :func:`genetic ancestry similarity estimation <pgscatalog.calc.AggregatedPGS.adjust>`.

    >>> from ._config import Config
    >>> score_path = Config.ROOT_DIR / "tests" / "data" / "aggregated_scores.txt.gz"
    >>> AggregatedPGS(path=score_path, target_name="hgdp")
    AggregatedPGS(path=PosixPath('.../aggregated_scores.txt.gz'))
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

    def adjust(self, *, ref_pc, target_pc, adjust_arguments=None):
        """Adjust a PGS based on genetic ancestry similarity estimations.

        :returns: :class:`AdjustResults`

        >>> from ._config import Config
        >>> from .principalcomponents import PrincipalComponents
        >>> related_path = Config.ROOT_DIR / "tests" / "data" / "ref.king.cutoff.id"
        >>> ref_pc = PrincipalComponents(pcs_path=[Config.ROOT_DIR / "tests" / "data" / "ref.pcs"], dataset="reference", psam_path=Config.ROOT_DIR / "tests" / "data" / "ref.psam", pop_type=PopulationType.REFERENCE, related_path=related_path)
        >>> target_pcs = PrincipalComponents(pcs_path=Config.ROOT_DIR / "tests" / "data" / "target.pcs", dataset="target", pop_type=PopulationType.TARGET)
        >>> score_path = Config.ROOT_DIR / "tests" / "data" / "aggregated_scores.txt.gz"
        >>> results = AggregatedPGS(path=score_path, target_name="hgdp").adjust(ref_pc=ref_pc, target_pc=target_pcs)
        >>> results.pgs.to_dict().keys()
        dict_keys(['SUM|PGS001229_hmPOS_GRCh38_SUM', 'SUM|PGS001229_hmPOS_GRCh38_AVG', 'percentile_MostSimilarPop|PGS001229_hmPOS_GRCh38_SUM', ...

        >>> results.models
        {'dist_empirical': {'PGS001229_hmPOS_GRCh38_SUM': {'EUR': {'percentiles': array([-1.04069000e+01, -7.94665080e+00, ...

        Write the adjusted results to a directory:

        >>> import tempfile, os
        >>> dout = tempfile.mkdtemp()
        >>> results.write(directory=dout)
        >>> sorted(os.listdir(dout))
        ['target_info.json.gz', 'target_pgs.txt.gz', 'target_popsimilarity.txt.gz']
        """
        if adjust_arguments is None:
            adjust_arguments = AdjustArguments()  # uses default values

        if ref_pc.pop_type != PopulationType.REFERENCE:
            raise ValueError("ref_pc argument has bad pop_type")

        if target_pc.pop_type != PopulationType.TARGET:
            raise ValueError("target_pc argument has bad pop_type")

        self._check_overlap(ref_pc=ref_pc, target_pc=target_pc)
        scorecols = list(self.df.columns)

        # join pgs + pca data
        target_df = target_pc.df.join(self.df.loc[self.target_name], on="IID")
        reference_df = ref_pc.df.join(self.df.loc["reference"], on="IID")

        assignment_threshold_p = choose_pval_threshold(adjust_arguments)

        if (ref_popcomp := ref_pc.npcs_popcomp) != (
            target_popcomp := target_pc.npcs_popcomp
        ):
            logger.warning(
                f"Reference PC: {ref_popcomp=} doesn't match {target_popcomp=}. Taking min()"
            )

        npcs_popcomp = min([ref_pc.npcs_popcomp, target_pc.npcs_popcomp])
        ancestry_ref, ancestry_target, compare_info = compare_ancestry(
            ref_df=reference_df,
            ref_pop_col=ref_pc.poplabel,
            ref_train_col="Unrelated",
            target_df=target_df,
            n_pcs=npcs_popcomp,
            method=adjust_arguments.method_compare,
            p_threshold=assignment_threshold_p,
        )

        reference_df = pd.concat([reference_df, ancestry_ref], axis=1)
        target_df = pd.concat([target_df, ancestry_target], axis=1)

        if (ref_norm := ref_pc.npcs_norm) != (target_norm := ref_pc.npcs_norm):
            logger.warning(
                f"Reference PC: {ref_norm=} doesn't match {target_norm=}. Taking min()"
            )

        npcs_norm = min([ref_pc.npcs_norm, target_pc.npcs_norm])

        # Adjust PGS values
        adjpgs_ref, adjpgs_target, adjpgs_models = pgs_adjust(
            reference_df,
            target_df,
            scorecols,
            ref_pc.poplabel,
            "MostSimilarPop",
            use_method=list(adjust_arguments.method_normalization),
            ref_train_col="Unrelated",
            n_pcs=npcs_norm,
        )
        adjpgs = pd.concat([adjpgs_ref, adjpgs_target], axis=0)

        reference_df["REFERENCE"] = True
        target_df["REFERENCE"] = False
        final_df = pd.concat([target_df, reference_df], axis=0)

        return AdjustResults(
            target_label=target_pc.dataset,
            models=adjpgs_models,
            model_meta=compare_info,
            pgs=adjpgs,
            pca=final_df,
            scorecols=scorecols,
        )


class PolygenicScore:
    """Represents the output of ``plink2 --score`` written to a file

    >>> from ._config import Config
    >>> import reprlib
    >>> score1 = Config.ROOT_DIR / "tests" / "data" / "cineca_22_additive_0.sscore.zst"
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
        """A generator that yields dataframe chunks."""
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
