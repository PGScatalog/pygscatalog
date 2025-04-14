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
                        index=["sampleset", "FID", "IID", "PGS"],
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
        # build sample IDs from (FID, IID)
        pca_ref_samples = set((x[1], x[2]) for x in ref_pc.df.index.tolist())
        pca_target_samples = set((x[1], x[2]) for x in target_pc.df.index.tolist())
        score_ref_samples = set(self.df.loc["reference"].index.tolist())
        score_target_samples = set(self.df.loc[self.target_name].index.tolist())

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
        >>> score_path = Config.ROOT_DIR / "tests" / "data" / "aggregated_scores.txt"
        >>> results = AggregatedPGS(path=score_path, target_name="hgdp").adjust(ref_pc=ref_pc, target_pc=target_pcs)
        >>> results.pgs.to_dict().keys()
        dict_keys(['SUM|PGS001229_hmPOS_GRCh38', 'percentile_MostSimilarPop|PGS001229_hmPOS_GRCh38', 'Z_MostSimilarPop|PGS001229_hmPOS_GRCh38', ...

        >>> results.models
        {'dist_empirical': {'PGS001229_hmPOS_GRCh38': {'EUR': {'percentiles': array([-1.04069000e+01, -7.94665080e+00, ...

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
        target_df = target_pc.df.join(self.df.loc[self.target_name], on=["FID", "IID"])
        reference_df = ref_pc.df.join(self.df.loc["reference"], on=["FID", "IID"])

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
    "{'DENOM': {('test', 'HG00096', 'HG00096'): 1564, ... 'PGS001229_22_SUM': {('test', 'HG00096', 'HG00096'): 0.54502, ...

    It's often helpful to combine PGS that were split per chromosome or by effect type:

    >>> aggregated_score = pgs1 + pgs2
    >>> aggregated_score  # doctest: +ELLIPSIS
    PolygenicScore(sampleset='test', path='(in-memory)')

    Once a score has been fully aggregated it can be helpful to recalculate an average:

    >>> aggregated_score.average()
    >>> aggregated_score.df  # doctest: +ELLIPSIS,+NORMALIZE_WHITESPACE
                                        PGS       SUM  DENOM       AVG
    sampleset FID     IID
    test      HG00096 HG00096  PGS001229_22  1.090040   3128  0.000348
              HG00097 HG00097  PGS001229_22  1.348802   3128  0.000431
    ...

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

    If a sampleset can't be inferred from argument or path, error:
    >>> PolygenicScore()
    Traceback (most recent call last):
    ...
    TypeError: Missing sampleset
    """

    def __init__(self, *, path=None, df=None, sampleset=None):
        try:
            self._path = pathlib.Path(path)
        except TypeError:
            self._path = None

        try:
            if sampleset is None:
                # infer missing sampleset from path
                self.sampleset = path.stem.split("_")[0]
            else:
                self.sampleset = sampleset
        except AttributeError:
            self.sampleset = sampleset
        finally:
            if self.sampleset is None:
                raise TypeError("Missing sampleset")

        self._df = df
        self._melted = False

    def __repr__(self):
        if self.path is None:
            path = repr("(in-memory)")
        else:
            path = repr(self.path)
        return f"{type(self).__name__}(sampleset={repr(self.sampleset)}, path={path})"

    def __add__(self, other):
        if isinstance(other, PolygenicScore):
            logger.info(f"Doing element-wise addition: {self} + {other}")
            sumdf = self.df.add(other.df, fill_value=0)
            return PolygenicScore(sampleset=self.sampleset, df=sumdf)
        else:
            return NotImplemented

    @property
    def df(self):
        if self._df is None:
            self._df = self.read()

        return self._df

    @property
    def path(self):
        return self._path

    def read(self):
        """Eagerly load a PGS into a pandas dataframe

        If the FID column can be missing from the input data:

        >>> from ._config import Config
        >>> from xopen import xopen
        >>> score1 = Config.ROOT_DIR / "tests" / "data" / "cineca_22_additive_0.sscore.zst"
        >>> with xopen(score1) as f:
        ...     f.readline().split()
        ['#IID', 'ALLELE_CT', 'DENOM', 'NAMED_ALLELE_DOSAGE_SUM', 'PGS001229_22_AVG', 'PGS001229_22_SUM']

        Then FID is set to IID:

        >>> PolygenicScore(sampleset="test", path=score1).read()  # doctest: +ELLIPSIS,+NORMALIZE_WHITESPACE
                                    DENOM  PGS001229_22_SUM
        sampleset FID     IID
        test      HG00096 HG00096   1564          0.545020
        ...
        """
        if self.path is None:
            raise ValueError("Missing path")

        df = pd.read_csv(
            self.path, sep="\t", converters={"FID": str, "IID": str}
        ).assign(sampleset=self.sampleset)

        if "#FID" not in df.columns:
            logger.warning("#FID column missing, setting FID == IID")
            # if FID is missing, IID starts with a hash
            df["#FID"] = df["#IID"]
            df = df.set_index(["sampleset", "#FID", "#IID"])
        elif all(df["#FID"] == 0):
            logger.warning("All FID column values missing (0), setting FID == IID")
            # FID column present, but missing data, IID doesn't start with a hash
            df["#FID"] = df["IID"]
            df = df.set_index(["sampleset", "#FID", "IID"])
        else:
            logger.info("#FID column detected")
            df = df.set_index(["sampleset", "#FID", "IID"])

        df.index.names = ["sampleset", "FID", "IID"]
        df = df[_select_agg_cols(df.columns)]
        return df

    def average(self):
        """Update the dataframe with a recalculated average."""
        logger.info("Recalculating average")
        if not self._melted:
            self.melt()

        df = self.df
        df["AVG"] = df.SUM / df.DENOM
        self._df = df

    def melt(self):
        """Update the dataframe with a melted version (wide format to long format)"""
        logger.info("Melting dataframe from wide to long format")
        df = self.df.melt(
            id_vars=["DENOM"],
            value_name="SUM",
            var_name="PGS",
            ignore_index=False,
        )
        # e.g. PGS000822_SUM -> PGS000822
        df["PGS"] = df["PGS"].str.replace("_SUM", "")
        # melted chunks need a consistent column order
        self._df = df[["PGS", "SUM", "DENOM"]]
        self._melted = True

    def write(self, outdir, split=False):
        """Write PGS to a compressed TSV"""
        outdir = pathlib.Path(outdir)

        if not self._melted:
            self.melt()

        df = self.df

        if split:
            logger.info("Writing results split by sampleset")

            for sampleset, group in df.groupby("sampleset"):
                fout = outdir / f"{sampleset}_pgs.txt.gz"
                df.to_csv(fout, sep="\t", compression="gzip", mode="a")
        else:
            logger.info("Writing combined results (aggregated_scores.txt.gz)")
            fout = outdir / "aggregated_scores.txt.gz"
            df.to_csv(fout, sep="\t", compression="gzip", mode="a")


def _select_agg_cols(cols):
    """Select aggregatable columns"""
    keep_cols = ["DENOM"]
    return [
        x
        for x in cols
        if (x.endswith("_SUM") and (x != "NAMED_ALLELE_DOSAGE_SUM")) or (x in keep_cols)
    ]
