import enum
import itertools
import logging

from .ancestry import read

logger = logging.getLogger(__name__)


class PopulationType(enum.Enum):
    TARGET = "target"
    REFERENCE = "reference"


class PrincipalComponents:
    """
    This class represents principal components data calculated by fraposa-pgsc

    >>> from ._config import Config
    >>> related_path = Config.ROOT_DIR / "tests" / "ref.king.cutoff.id"
    >>> psam_path = Config.ROOT_DIR / "tests" / "ref.psam"
    >>> ref_pc = PrincipalComponents(pcs_path=[Config.ROOT_DIR / "tests" / "ref.pcs"], dataset="reference", psam_path=psam_path, related_path=related_path, pop_type=PopulationType.REFERENCE)
    >>> ref_pc
    PrincipalComponents(dataset='reference', pop_type=PopulationType.REFERENCE, pcs_path=[PosixPath('.../pgscatalog.calclib/tests/ref.pcs')], psam_path=PosixPath('.../pgscatalog.calclib/tests/ref.psam'))
    >>> ref_pc.df.to_dict()
    {'PC1': {('reference', 'HG00096'): -23.8212, ('reference', 'HG00097'): -24.8106, ...

    >>> target_pcs = PrincipalComponents(pcs_path=Config.ROOT_DIR / "tests" / "target.pcs", dataset="target", pop_type=PopulationType.TARGET)
    >>> target_pcs
    PrincipalComponents(dataset='target', pop_type=PopulationType.TARGET, pcs_path=[PosixPath('.../pgscatalog.calclib/tests/target.pcs')], psam_path=None)
    >>> target_pcs.df.to_dict()
    {'PC1': {('target', 'HGDP00001'): -18.5135, ('target', 'HGDP00003'): -18.8314, ...
    """

    def __init__(
        self, pcs_path, dataset, pop_type, psam_path=None, related_path=None, **kwargs
    ):
        try:
            self.pcs_path = list(itertools.chain(pcs_path))  # handle a list of paths
        except TypeError:
            self.pcs_path = [pcs_path]  # or a single path

        self.dataset = dataset
        self._psam_path = psam_path
        self._related_path = related_path
        self._pop_type = pop_type

        if self.pop_type == PopulationType.REFERENCE:
            if self.psam_path is None or self.related_path is None:
                raise ValueError("Reference requires psam_path and related_path")

        self._npcs_popcomp = kwargs.get("npcs_popcomp", 5)
        self._npcs_norm = kwargs.get("npcs_normalization", 4)

        self._df = None

        if self.pop_type == PopulationType.REFERENCE:
            # Population labels in REFERENCE psam to use for assignment
            self._poplabel = kwargs.get("ref_label", "SuperPop")
        else:
            self._poplabel = None

    def __repr__(self):
        return f"PrincipalComponents(dataset={self.dataset!r}, pop_type={self.pop_type}, pcs_path={self.pcs_path!r}, psam_path={self.psam_path!r})"

    @property
    def pop_type(self):
        return self._pop_type

    @property
    def psam_path(self):
        return self._psam_path

    @property
    def related_path(self):
        return self._related_path

    @property
    def poplabel(self):
        return self._poplabel

    @property
    def max_pcs(self):
        """The maximum number of PCs used in calculations"""
        return max([10, self._npcs_popcomp, self._npcs_norm])

    @property
    def npcs_popcomp(self):
        """Number of PCs used for population comparison (default = 5)"""
        return self._npcs_popcomp

    @npcs_popcomp.setter
    def npcs_popcomp(self, value):
        if 1 <= value <= 20:
            self._npcs_popcomp = int(value)
        else:
            raise ValueError("Must be integer between 1 and 20")

    @property
    def npcs_norm(self):
        """Number of PCs used for population NORMALIZATION (default = 4)"""
        return self._npcs_norm

    @npcs_norm.setter
    def npcs_norm(self, value):
        if 1 <= value <= 20:
            self._npcs_popcomp = int(value)
        else:
            raise ValueError("Must be integer between 1 and 20")

    @property
    def df(self):
        if self._df is None:
            df = read.read_pcs(
                loc_pcs=self.pcs_path,
                dataset=self.dataset,
                loc_related_ids=self.related_path,
                nPCs=self.max_pcs,
            )
            if self.pop_type == PopulationType.REFERENCE:
                df = read.extract_ref_psam_cols(
                    loc_psam=self.psam_path,
                    dataset=self.dataset,
                    df_target=df,
                    keepcols=self._poplabel,
                )
            self._df = df

        if self._df.shape[0] < 100 and self.pop_type == PopulationType.REFERENCE:
            logger.critical(
                "Error: too few reference panel samples. This is an arbitrary threshold "
                "for input QC; however, it is inadvisable to run this analysis with limited "
                "reference panel diversity as empirical percentiles are calculated."
            )
            raise ValueError

        return self._df
