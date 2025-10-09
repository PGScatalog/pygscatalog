from __future__ import annotations

import logging
import pathlib
import shutil
import sqlite3
from abc import ABC, abstractmethod
from functools import cached_property
from typing import TYPE_CHECKING

import pysam

from ._bgen import bgen_buffer_variants, bgen_get_sample_list
from ._vcf import vcf_buffer_variants, vcf_get_sample_list
from .genomefiletypes import GenomeFileType

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    from .targetvariant import TargetVariant
    from .types import Pathish


class GenomeFileHandler(ABC):
    def __init__(self, path: Pathish, cache_dir: Pathish):
        self._target_path = pathlib.Path(path).resolve()
        self._cache_dir = pathlib.Path(cache_dir).resolve()

        if not self._target_path.exists():
            raise FileNotFoundError(self._target_path)

        if not self._cache_dir.is_dir():
            raise NotADirectoryError(self._cache_dir)

    @property
    @abstractmethod
    def file_type(self) -> GenomeFileType: ...

    @property
    @abstractmethod
    def chroms(self) -> list[str]: ...

    @property
    @abstractmethod
    def samples(self) -> list[str]: ...

    @property
    def target_path(self) -> pathlib.Path:
        return self._target_path

    @property
    @abstractmethod
    def index_path(self) -> pathlib.Path: ...

    @index_path.setter
    @abstractmethod
    def index_path(self, path: pathlib.Path) -> None: ...

    @abstractmethod
    def query_variants(
        self, positions: Sequence[tuple[str, int]]
    ) -> Iterable[TargetVariant]: ...


class VCFHandler(GenomeFileHandler):
    def __init__(self, path: Pathish, cache_dir: Pathish):
        if not pathlib.Path(path).suffixes[-2:] == [".vcf", ".gz"]:
            raise ValueError(f"{self.target_path} is not a .vcf.gz file")

        super().__init__(path=path, cache_dir=cache_dir)

        csi_path = pathlib.Path(str(self._target_path) + ".csi")
        if csi_path.exists():
            self._index_path = csi_path
        else:
            tbi_path = pathlib.Path(str(self._target_path) + ".tbi")
            if tbi_path.exists():
                self._index_path = tbi_path
            else:
                raise FileNotFoundError(f"{csi_path} or {tbi_path}")

    @property
    def index_path(self) -> pathlib.Path:
        return self._index_path

    @index_path.setter
    def index_path(self, path: pathlib.Path) -> None:
        if not path.resolve().exists():
            raise FileNotFoundError(path)
        self._index_path = path

    @property
    def file_type(self) -> GenomeFileType:
        return GenomeFileType.VCF

    @property
    def chroms(self) -> list[str]:
        with pysam.VariantFile(
            str(self.target_path), mode="r", drop_samples=True
        ) as vcf:
            return [str(x) for x in vcf.header.contigs]

    def query_variants(
        self,
        positions: Sequence[tuple[str, int]],
    ) -> Iterable[TargetVariant]:
        return vcf_buffer_variants(
            position_batch=positions,
            cache_dir=self._cache_dir,
            target_path=self._target_path,
        )

    @property
    def samples(self) -> list[str]:
        return vcf_get_sample_list(target_path=self.target_path)


class BgenFileHandler(GenomeFileHandler):
    def __init__(self, path: Pathish, cache_dir: Pathish, sample_file: Pathish):
        super().__init__(path=path, cache_dir=cache_dir)
        self._sample_file = pathlib.Path(sample_file)
        self._index_path = self._target_path.with_suffix(".bgen.bgi")

        if not self._sample_file.exists():
            raise FileNotFoundError(self._sample_file)

        if not self._index_path.exists():
            raise FileNotFoundError(self._index_path)

    @property
    def sample_file(self) -> pathlib.Path:
        return self._sample_file

    @property
    def file_type(self) -> GenomeFileType:
        return GenomeFileType.BGEN

    @property
    def index_path(self) -> pathlib.Path:
        index_path = self.target_path.with_name(
            f"{self.target_path.name}.bgi"
        ).resolve()
        if not index_path.exists():
            raise FileNotFoundError(f"{index_path} does not exist")
        return index_path

    @index_path.setter
    def index_path(self, path: pathlib.Path) -> None:
        if not path.resolve().exists():
            raise FileNotFoundError(path)
        self._index_path = path

    @cached_property
    def chroms(self) -> list[str]:
        with sqlite3.connect(self._index_path) as conn:
            rows = conn.execute("SELECT DISTINCT chromosome FROM Variant").fetchall()
        return [row[0] for row in rows]

    def query_variants(
        self,
        positions: Sequence[tuple[str, int]],
    ) -> Iterable[TargetVariant]:
        # important to copy the bgen index to the local cache
        # indexes are often mounted read-only and they need to be modified for queries
        # (they're just sqlite databases! very cool)
        self.copy_index_to_cache()
        return bgen_buffer_variants(
            cache_dir=self._cache_dir,
            position_batch=positions,
            target_path=self.target_path,
            idx_path=self._index_path,
            target_chroms=self.chroms,
        )

    def copy_index_to_cache(self) -> None:
        """Copy the index file to the cache directory.

        Helpful when genomes are mounted read only. Only used for bgen files.
        """
        cache_index_path = self._cache_dir / "tmp" / self.index_path.name
        cache_index_path.parent.mkdir(parents=True, exist_ok=True)
        # VCF indexes aren't staged to the cache
        logger.info(f"Staging bgi index {self.index_path} to {cache_index_path}")
        if not cache_index_path.exists():
            logger.info("Copying bgi index to cache")
            # copyfile prevents preserving file attributes (read-only)
            self._index_path = shutil.copyfile(self._index_path, cache_index_path)
        else:
            logger.info("bgi index already in cache")
            self._index_path = cache_index_path.resolve()

    @property
    def samples(self) -> list[str]:
        return bgen_get_sample_list(
            target_path=self.target_path, sample_path=self.sample_file
        )


def get_file_handler(
    path: Pathish,
    *,
    cache_dir: Pathish,
    sample_file: Pathish | None = None,
) -> GenomeFileHandler:
    path = pathlib.Path(path)

    if path.suffixes[-2:] == [".vcf", ".gz"]:
        return VCFHandler(path=path, cache_dir=cache_dir)
    if path.suffix == ".bgen":
        if sample_file is None:
            raise ValueError("BGEN files require a sample file")
        return BgenFileHandler(path=path, cache_dir=cache_dir, sample_file=sample_file)
    raise ValueError(f"Unsupported genome file type: {path}")
