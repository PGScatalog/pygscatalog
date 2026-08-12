"""MATCH: match scoring file variants against one chromosome's shards.

Staged with every shard of one chromosome, so it matches against a complete variant
list and no shard has a partial view. It reads the cache's 1D variant metadata
only, never the 3D genotype arrays, and it reads the sampleset name from the cache
rather than being told it.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, cast

import polars as pl

from pgscatalog.calc.lib.cache.genotypecache import GenotypeCache, check_cache_shape

from . import stage
from ._matchvariants import update_match_table
from .scoredb import ScoreDatabase

if TYPE_CHECKING:
    import duckdb
    import zarr

    from pgscatalog.calc.lib.types import Pathish, PathishList


logger = logging.getLogger(__name__)

# columns in the targetvariants table
_TARGET_VARIANTS_COLUMNS = [
    "geno_index",
    "chr_name",
    "chr_pos",
    "ref",
    "alts",
    "variant_id",
    "zarr_group",
    "sampleset",
]


class MatchStage:
    """One MATCH task: the shards of one chromosome, and the database it matches in.

    Parameters
    ----------
    out_dir : Pathish
        The task's output directory. `matches.parquet` and `samples.parquet` are
        published at its root; the cache and the database live in its scratch.
    max_memory_gb : float
        Maximum RAM duckdb may use before intermediate results spill to disk.
    threads : int
        Maximum number of threads duckdb may use.
    """

    def __init__(self, out_dir: Pathish, max_memory_gb: float, threads: int):
        scratch = stage.scratch_dir(out_dir)

        # the cache owns the zarr hierarchy, so this stage asks it what it holds
        # instead of deriving paths. MATCH is staged with a cache, so it reads
        # the sampleset name from it
        self._cache = GenotypeCache(path=scratch / "genotypes.zarr")

        # every query this stage runs goes through the database, so the resource
        # limits the CLI was given reach all of them
        self._db = ScoreDatabase(
            path=scratch / "pgs.db",
            max_memory_gb=max_memory_gb,
            threads=threads,
        )

    def load_scores(self, scorefile_paths: Pathish | PathishList) -> None:
        """Copy scores processed with pgscatalog-format into a duckDB database table
        'scorevariants'
        """
        stage.load_scores(db=self._db, scorefile_paths=scorefile_paths)

    def _load_variants_from_zarr(self, conn: duckdb.DuckDBPyConnection) -> None:
        """Load variant metadata from zarr groups into targetvariants table

        Variant metadata is stored as 1D arrays in "meta" groups

        The cache owns the hierarchy, so it hands over each shard's metadata; this
        builds a dataframe per shard and inserts it
        """
        root: zarr.Group = self._cache.root
        sampleset: str = self._cache.sampleset

        conn.sql("""
        CREATE TABLE IF NOT EXISTS targetvariants (
            geno_index UINTEGER,
            chr_name TEXT,
            chr_pos BIGINT,
            ref TEXT,
            alts TEXT[],
            variant_id TEXT,
            zarr_group TEXT,
            sampleset TEXT
        );
        """)

        # iterate over each shard's metadata lazily to save memory, because these 1D
        # arrays can get big, especially when loaded into a polars DF
        for zarr_group, variants in self._cache.variant_metadata():
            logger.info(f"Inserting variants from zarr {zarr_group=}")
            _df = (
                (
                    pl.DataFrame(variants)
                    # replace all empty strings with nulls
                    .with_columns(pl.col(pl.String).replace("", None))
                )
                # re-split the alts field
                .with_columns(
                    [
                        pl.lit(zarr_group).alias("zarr_group"),
                        pl.lit(sampleset).alias("sampleset"),
                        pl.col("alts").str.split("/"),
                    ]
                )
                # row index is really important!
                # corresponds to row numbers in associated genotype array (0-based)
                .with_row_index("geno_index")
                # correct order for ingest
                .select(_TARGET_VARIANTS_COLUMNS)
            )

            # check the number of variants in metadata matches the genotype array shape
            sub_group = cast("zarr.Group", root[zarr_group])
            n_geno_variants: int = cast("zarr.Array", sub_group["genotypes"]).shape[0]

            if n_geno_variants == (n_meta_variants := _df.shape[0]):
                logger.info(
                    "Number of variants in metadata matches number of variants in "
                    "genotypes array"
                )
            else:
                logger.critical(f"Loading {zarr_group=} {sampleset=} failed")
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
                    zarr_group,
                    sampleset
                FROM _df""")
            logger.info(f"{_df.shape=} finished loading into targetvariants table")

        conn.sql("""
        -- data integrity check after loading is finished
        -- zarr_group is the full path to the shard, so it implies the sampleset
        ALTER TABLE targetvariants ADD PRIMARY KEY (variant_id, zarr_group);
        """)

    def match_variants(
        self,
        match_ambiguous: bool = False,
        match_multiallelic: bool = False,
        match_flipped: bool = True,
    ) -> str:
        """Match scoring file variants against one chromosome's shards.

        No `min_overlap` here: the match rate is genome-wide, so only
        MATCH_COMBINE can decide it.

        Returns the chromosome that was matched.
        """
        logger.info("Matching scoring file variants in target genomes")

        with self._db.connect() as conn:
            # read variant information from zarr arrays into DB
            self._load_variants_from_zarr(conn)

            # check the shards can be scored independently before doing the work
            chromosome = check_cache_shape(
                conn=conn,
                shards=self._cache.shards,
                sampleset=self._cache.sampleset,
            )

            # match score variants to zarr variants
            update_match_table(
                conn=conn,
                match_ambiguous=match_ambiguous,
                match_multiallelic=match_multiallelic,
                match_flipped=match_flipped,
                sampleset=self._cache.sampleset,
            )
            return chromosome

    def export_match_table(self, out_path: Pathish) -> None:
        """Write this chromosome's match table for MATCH_COMBINE to pool"""
        logger.info(f"Exporting match table to {str(out_path)}")
        with self._db.connect() as conn:
            conn.execute(
                "COPY allele_match_table TO ? (FORMAT parquet, COMPRESSION zstd)",
                [str(out_path)],
            )

    def export_samples(self, out_path: Pathish, chromosome: str) -> None:
        """Write the chromosome's sample IDs so MATCH_COMBINE can check the union"""
        logger.info(f"Exporting samples to {str(out_path)}")
        samples = sorted(
            sample for shard in self._cache.shards for sample in shard.samples
        )
        with self._db.connect() as conn:
            conn.execute(
                """
                COPY (
                    SELECT $sampleset AS sampleset,
                           $chromosome AS chromosome,
                           UNNEST($samples) AS sample_id
                ) TO $path (FORMAT parquet, COMPRESSION zstd)
                """,
                {
                    "sampleset": self._cache.sampleset,
                    "chromosome": chromosome,
                    "samples": samples,
                    "path": str(out_path),
                },
            )
