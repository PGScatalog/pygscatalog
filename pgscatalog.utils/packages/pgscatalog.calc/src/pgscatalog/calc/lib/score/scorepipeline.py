from __future__ import annotations

import logging
import pathlib
from typing import TYPE_CHECKING, cast

import dask.config
import duckdb
import polars as pl
import zarr
import zarr.storage

from pgscatalog.calc.lib.scorefile import load_scoring_files

from ._dosage import (
    get_sample_ids,
    make_dosage_arrays,
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

if TYPE_CHECKING:
    from collections.abc import Generator

    import numpy as np
    import numpy.typing as npt

    from pgscatalog.calc.lib.types import Pathish, PathishList


logger = logging.getLogger(__name__)


class ScorePipeline:
    def __init__(
        self,
        max_memory_gb: float,
        threads: int,
        out_dir: Pathish,
        minimum_samples_for_impute: int = 50,
    ):
        self._out_dir = pathlib.Path(out_dir).resolve()

        # the output directory contains results, so don't overwrite stuff
        if not self._out_dir.is_dir():
            raise NotADirectoryError(f"{str(self._out_dir)} must be a directory")

        self._n_minimum_samples_for_impute = minimum_samples_for_impute

        # note: duckdb will use 90% of available memory and all threads by default
        self._max_memory_gb = f"{max_memory_gb}GB"
        self._threads = threads

    @property
    def out_dir(self) -> pathlib.Path:
        return self._out_dir

    @property
    def db_path(self) -> pathlib.Path:
        return self._out_dir / "pgs.db"

    @property
    def genotypes_root_path(self) -> str:
        return str((self._out_dir / "genotypes.zarr").resolve())

    @property
    def sampleset(self) -> str:
        """Read sampleset name from the zarr directory"""
        zarr = pathlib.Path(self.genotypes_root_path)
        samplesets = [p for p in zarr.iterdir() if p.is_dir()]
        if len(samplesets) > 1:
            raise NotImplementedError(
                f"One than one sampleset not supported, found {samplesets}"
            )
        if len(samplesets) == 0:
            raise FileNotFoundError(
                "zarr archive looks badly structured, can't find sampleset group"
            )
        return str(samplesets[0].name)

    @property
    def target_variants_columns(self) -> list[str]:
        # columns in the targetvariants table
        return [
            "geno_index",
            "chr_name",
            "chr_pos",
            "ref",
            "alts",
            "variant_id",
            "filename",
            "sampleset",
        ]

    def load_scores(self, scorefile_paths: Pathish | PathishList) -> None:
        """Copy scores processed with pgscatalog-format into a duckDB database table
        'scorevariants'
        """
        if not isinstance(scorefile_paths, list):
            scorefile_paths = [scorefile_paths]

        load_scoring_files(
            db_path=self.db_path,
            scorefile_paths=scorefile_paths,
            max_memory_gb=self._max_memory_gb,
            threads=self._threads,
        )

    def _load_variants_from_zarr(self, conn: duckdb.DuckDBPyConnection) -> None:
        """Load variant metadata from zarr groups into targetvariants table

        Variant metadata is stored as 1D arrays in "meta" groups

        This function traverses the hierarchy, builds a dataframe, and inserts it
        """
        store = zarr.storage.LocalStore(root=self.genotypes_root_path)
        # load the sampleset group ready for traversing
        group: zarr.Group = zarr.open_group(store=store, path=self.sampleset, mode="r")

        conn.sql("""
        CREATE TABLE IF NOT EXISTS targetvariants (
            geno_index UINTEGER,
            chr_name TEXT,
            chr_pos BIGINT,
            ref TEXT,
            alts TEXT[],
            variant_id TEXT,
            filename TEXT,
            sampleset TEXT
        );
        """)

        # when traversing the group, iterate over each meta group lazily to save memory
        # because these 1D arrays can get big, especially when loaded into a polars DF
        for filename, variants in get_variants_from_zarr(group):
            logger.info(f"Inserting variants from zarr {filename=}")
            _df = (
                (
                    pl.DataFrame(variants)
                    # replace all empty strings with nulls
                    .with_columns(pl.col(pl.String).replace("", None))
                )
                # re-split the alts field
                .with_columns(
                    [
                        pl.lit(filename).alias("filename"),
                        pl.lit(self.sampleset).alias("sampleset"),
                        pl.col("alts").str.split("/"),
                    ]
                )
                # row index is really important!
                # corresponds to row numbers in associated genotype array (0-based)
                .with_row_index("geno_index")
                # correct order for ingest
                .select(self.target_variants_columns)
            )

            # check the number of variants in metadata matches the genotype array shape
            sub_group = cast("zarr.Group", group[filename])
            n_geno_variants: int = cast("zarr.Array", sub_group["genotypes"]).shape[0]

            if n_geno_variants == (n_meta_variants := _df.shape[0]):
                logger.info(
                    "Number of variants in metadata matches number of variants in "
                    "genotypes array"
                )
            else:
                logger.critical(f"Loading {filename=} {self.sampleset=} failed")
                raise ValueError(f"{n_geno_variants=} but {n_meta_variants=}")

            # it's OK to insert the variant data to the database now
            conn.sql("""
                INSERT INTO targetvariants
                SELECT
                    geno_index,
                    chr_name,
                    chr_pos,
                    ref,
                    alts,
                    variant_id,
                    filename,
                    sampleset
                FROM _df""")
            logger.info(f"{_df.shape=} finished loading into targetvariants table")

        conn.sql("""
        -- data integrity check after loading is finished
        ALTER TABLE targetvariants ADD PRIMARY KEY (variant_id, filename, sampleset);
        """)

    def match_variants(
        self,
        match_ambiguous: bool = False,
        match_multiallelic: bool = False,
        min_overlap: float = 0.75,
    ) -> None:
        logger.info("Matching scoring file variants in target genomes")

        if not (0 < min_overlap <= 1):
            raise ValueError(
                f"min_overlap must be positive and less than 1, got {min_overlap}"
            )

        with duckdb.connect(
            self.db_path,
            config={
                "max_memory": self._max_memory_gb,
                "threads": str(self._threads),
            },
        ) as conn:
            # read variant information from zarr arrays into DB
            self._load_variants_from_zarr(conn)
            # match score variants to zarr variants
            update_match_table(
                conn=conn,
                match_ambiguous=match_ambiguous,
                match_multiallelic=match_multiallelic,
                sampleset=self.sampleset,
            )
            # create match log
            add_match_log(
                conn=conn,
                min_overlap=min_overlap,
                sampleset=self.sampleset,
            )

    def export_full_match_log(self, out_directory: Pathish) -> None:
        logger.info(f"Exporting full match log to {str(out_directory)}")
        pathlib.Path(out_directory).mkdir(parents=True, exist_ok=True)

        with duckdb.connect(
            self.db_path,
            config={"max_memory": self._max_memory_gb, "threads": str(self._threads)},
        ) as conn:
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
            config={"max_memory": self._max_memory_gb, "threads": str(self._threads)},
        ) as conn:
            (
                conn.table("summary_log_table")
                .order("sampleset,accession,match_summary")
                .write_csv(file_name=str(out_path))
            )

        logger.info("Summary match log exported")

    def export_scores(self, out_path: Pathish) -> None:
        write_scores(
            db_path=self.db_path,
            out_dir=out_path,
            max_memory_gb=self._max_memory_gb,
            threads=self._threads,
        )

    @property
    def genotypes_root(self) -> zarr.Group:
        store = zarr.storage.LocalStore(self._out_dir / "genotypes.zarr")
        return zarr.open_group(store=store, mode="r")

    @property
    def samplesets(self) -> list[str]:
        """Get samplesets in the zarr store"""
        samplesets = list(self.genotypes_root.keys())
        if len(samplesets) > 1:
            raise NotImplementedError
        return samplesets

    def calculate_scores(self) -> None:
        create_score_table(db_path=self.db_path)

        for sampleset in self.samplesets:
            pgs_store = zarr.storage.LocalStore(self.out_dir)
            pgs_group = zarr.open_group(
                store=pgs_store, mode="w", path=f"pgs/{sampleset}"
            )

            # get accessions to include in score calculation
            # (based on variant match thresholds)
            accessions_to_calculate = get_ok_accessions(
                db_path=self.db_path, sampleset=sampleset
            )

            # get wide effect weights (a matrix) for each zarr genotype group
            grouped_dfs: dict[str, pl.DataFrame] = store_group_weight_arrays(
                db_path=self.db_path,
                sampleset=sampleset,
                pgs_group=pgs_group,
                accessions=accessions_to_calculate,
            )

            # calculate and adjust dosage with dask
            is_missing_zarr, dosage_zarr = make_dosage_arrays(
                genotypes_group=self.genotypes_root,
                df_groups=grouped_dfs,
                pgs_store=pgs_store,
                sampleset=sampleset,
                db_path=self.db_path,
                n_workers=self._threads,
                n_minimum_impute=self._n_minimum_samples_for_impute,
            )

            with dask.config.set(scheduler="threads", num_workers=self._threads):
                # calculate weighted sum of scores
                effect_weights: zarr.Array = cast(
                    "zarr.Array", pgs_group["weight_matrix"]
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
                    sample_ids=get_sample_ids(
                        group=self.genotypes_root, sampleset=sampleset
                    ),
                    is_missing_array=is_missing_zarr,
                )

            # insert scores into duckdb for final cleanup and export
            insert_scores(
                db_path=self.db_path,
                _score_df=_score_df,
            )
            logger.info(f"Score calculation {sampleset=} finished")


def get_variants_from_zarr(
    group: zarr.Group,
) -> Generator[tuple[str, dict[str, npt.NDArrayLike]], None, None]:
    # traverse a zarr group and extract all 1D meta arrays
    # returning a dataframe structure
    for name, subgrp in group.groups():
        if name == "meta":
            meta = {}
            for arr_name, arr in subgrp.arrays():
                # loads the numpy array into memory now
                meta[arr_name] = arr[:]
            # meta is ready for conversion to a DF
            # (key = column name, value = 1D numpy vector/rows)
            yield group.basename, meta
        else:
            # recurse the group tree
            yield from get_variants_from_zarr(subgrp)
