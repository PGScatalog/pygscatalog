from __future__ import annotations

import logging
import pathlib
from collections.abc import Sequence
from typing import TYPE_CHECKING, cast

import dask.config
import duckdb
import numpy as np
import polars as pl
import zarr

from ._dosage import (
    store_dosage_from_chunks,
)
from ._matchlog import add_match_log, get_ok_accessions
from ._matchvariants import update_match_table
from ._pgs import (
    calculate_score_statistics,
    calculate_scores,
    create_score_table,
    insert_scores,
    write_scores,
)
from ._weight_matrix import store_group_weight_arrays
from .scorefile import load_scoring_files

if TYPE_CHECKING:
    import numpy.typing as npt

    from .types import Pathish, PathishList


logger = logging.getLogger(__name__)


# # TODO: add state checks
# class ScoreDBState(enum.Enum):
#     INIT = enum.auto()
#     SCORES_LOADED = enum.auto()
#     VARIANTS_MATCHED = enum.auto()
#     SCORES_CALCULATED = enum.auto()


class ScorePipeline:
    def __init__(
        self,
        cache_dir: Pathish,
        max_memory_gb: float,
        threads: int,
        db_path: Pathish | None = None,
        minimum_samples_for_impute: int = 50,
    ):
        self._cache_dir = pathlib.Path(cache_dir)

        if db_path is None:
            self._db_path = self._cache_dir / "scores.db"
        else:
            self._db_path = pathlib.Path(db_path)

        if self._db_path.exists():
            logger.warning("Score database file exists, overwriting")
            self._db_path.unlink()

        self._n_minimum_samples_for_impute = minimum_samples_for_impute
        # note: duckdb will use 90% of available memory and all threads by default
        self._max_memory_gb = f"{max_memory_gb}GB"
        self._threads = threads

    @property
    def db_path(self) -> pathlib.Path:
        return self._db_path

    def load_scores(self, scorefile_paths: Pathish | PathishList) -> None:
        """Copy scores processed with pgscatalog-format into a duckDB database table
        'scorevariants'
        """
        if not isinstance(scorefile_paths, list):
            scorefile_paths = [scorefile_paths]

        load_scoring_files(
            db_path=self._db_path,
            scorefile_paths=scorefile_paths,
            max_memory_gb=self._max_memory_gb,
            threads=self._threads,
        )

    def _attach_target_variants(self, conn: duckdb.DuckDBPyConnection) -> None:
        """Attach the target variants database (in the cache directory)

        Attachments are scoped per-sesssion, so this has to be called when working with
        target variants
        """
        variants_db_path = str((self._cache_dir / "variants.db").resolve())
        conn.sql(
            f"ATTACH IF NOT EXISTS '{variants_db_path}' AS target_variants (READ_ONLY)"
        )

    def match_variants(
        self,
        match_ambiguous: bool = False,
        match_multiallelic: bool = False,
        min_overlap: float = 0.75,
    ) -> None:
        logger.info("Matching scoring file variants in target genomes")

        if 0 >= min_overlap > 1:
            raise ValueError(
                f"min_overlap must be positive and less than 1, got {min_overlap}"
            )

        for sampleset in self.samplesets:
            with duckdb.connect(
                self.db_path,
                config={"max_memory": self._max_memory_gb, "threads": self._threads},
            ) as conn:
                self._attach_target_variants(conn)
                update_match_table(
                    conn=conn,
                    match_ambiguous=match_ambiguous,
                    match_multiallelic=match_multiallelic,
                    sampleset=sampleset,
                )

                add_match_log(
                    conn=conn,
                    min_overlap=min_overlap,
                    sampleset=sampleset,
                )

    def export_full_match_log(self, out_directory: Pathish) -> None:
        logger.info(f"Exporting full match log to {str(out_directory)}")
        pathlib.Path(out_directory).mkdir(parents=True, exist_ok=True)

        with duckdb.connect(
            self.db_path,
            config={"max_memory": self._max_memory_gb, "threads": self._threads},
        ) as conn:
            self._attach_target_variants(conn)
            conn.table("score_log_table").order("accession, row_nr").write_csv(
                file_name=str(out_directory),
                partition_by=["sampleset", "accession", "chr_name"],
                write_partition_columns=True,
                compression="gzip",
                overwrite=True,
            )
        logger.info("Full match log exported")

    def export_summary_match_log(self, out_path: Pathish) -> None:
        logger.info(f"Exporting summary match log to {str(out_path)}")
        pathlib.Path(out_path).parent.mkdir(parents=True, exist_ok=True)

        with duckdb.connect(
            self.db_path,
            config={"max_memory": self._max_memory_gb, "threads": self._threads},
        ) as conn:
            conn.table("summary_log_table").write_csv(file_name=str(out_path))

        logger.info("Summary match log exported")

    def export_scores(self, out_path: Pathish) -> None:
        write_scores(
            db_path=self.db_path,
            out_dir=out_path,
            max_memory_gb=self._max_memory_gb,
            threads=self._threads,
        )

    @property
    def samplesets(self) -> list[str]:
        """Get samplesets in the zarr store"""
        store = zarr.storage.LocalStore(self._cache_dir / "gts")
        root = zarr.open_group(store=store, mode="r")
        samplesets = list(root.keys())
        if len(samplesets) > 1:
            raise NotImplementedError
        else:
            return samplesets

    def calculate_scores(self) -> None:
        create_score_table(db_path=self.db_path)

        for sampleset in self.samplesets:
            store = zarr.storage.LocalStore(self._cache_dir)
            pgs_group = zarr.open_group(store=store, mode="w", path=f"pgs/{sampleset}")

            accessions_to_calculate = get_ok_accessions(
                db_path=self.db_path, sampleset=sampleset
            )

            # metadata df for each zarr group
            grouped_dfs: dict[str, pl.DataFrame] = store_group_weight_arrays(
                db_path=self.db_path,
                sampleset=sampleset,
                zarr_group=pgs_group,
                accessions=accessions_to_calculate,
            )

            # store dosage in zarr and adjust it
            is_missing_zarr, dosage_zarr = store_dosage_from_chunks(
                df_groups=grouped_dfs,
                store=store,
                sampleset=sampleset,
                db_path=self.db_path,
                n_workers=self._threads,
                n_minimum_impute=self._n_minimum_samples_for_impute,
            )

            with dask.config.set(scheduler="threads", num_workers=self._threads):
                effect_weights: zarr.Array = cast(
                    zarr.Array, pgs_group["weight_matrix"]
                )
                scores: npt.NDArray[np.float64] = calculate_scores(
                    dosage_array=dosage_zarr, effect_weights=effect_weights
                ).compute()

                _score_df = calculate_score_statistics(
                    db_path=self.db_path,
                    accessions=accessions_to_calculate,
                    dosage_array=dosage_zarr,
                    score=scores,
                    sampleset=sampleset,
                    sample_ids=get_sample_ids(store=store, sampleset=sampleset),
                    is_missing_array=is_missing_zarr,
                )

            insert_scores(
                db_path=self.db_path,
                _score_df=_score_df,
            )
            logger.info(f"Score calculation {sampleset=} finished")


def get_sample_ids(store: zarr.storage.StoreLike, sampleset: str) -> list[str]:
    group = zarr.open_group(store=store, path=f"gts/{sampleset}", mode="r")
    # it's definitely a list of strings
    samples = cast(Sequence[str], group.attrs["samples"])
    return [str(x) for x in samples]
