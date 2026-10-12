from __future__ import annotations

import itertools
import logging
import pathlib
from typing import TYPE_CHECKING, Any, cast

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
from ._matchvariants import (
    _setup_enums,
    create_allele_match_table,
    update_match_table,
)
from ._pgs import (
    calculate_score_statistics,
    calculate_scores,
    create_score_table,
    insert_scores,
    write_scores,
)
from ._weight_matrix import store_group_weight_arrays

if TYPE_CHECKING:
    from collections.abc import Generator, Sequence

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
        sampleset: str | None = None,
    ):
        # MATCH_COMBINE and SCORE_COMBINE have no genotypes.zarr to read the name
        # from, so they pass it in
        self._sampleset = sampleset
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
        if self._sampleset is not None:
            return self._sampleset

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
            "zarr_group",
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
        root: zarr.Group = zarr.open_group(store=store, mode="r")
        # load the sampleset group ready for traversing
        group: zarr.Group = cast("zarr.Group", root[self.sampleset])

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

        # when traversing the group, iterate over each meta group lazily to save memory
        # because these 1D arrays can get big, especially when loaded into a polars DF
        for zarr_group, variants in get_variants_from_zarr(group):
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
            sub_group = cast("zarr.Group", root[zarr_group])
            n_geno_variants: int = cast("zarr.Array", sub_group["genotypes"]).shape[0]

            if n_geno_variants == (n_meta_variants := _df.shape[0]):
                logger.info(
                    "Number of variants in metadata matches number of variants in "
                    "genotypes array"
                )
            else:
                logger.critical(f"Loading {zarr_group=} {self.sampleset=} failed")
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

    def _connect(self) -> duckdb.DuckDBPyConnection:
        return duckdb.connect(
            self.db_path,
            config={
                "max_memory": self._max_memory_gb,
                "threads": str(self._threads),
            },
        )

    def match_variants(
        self,
        match_ambiguous: bool = False,
        match_multiallelic: bool = False,
    ) -> str:
        """Match scoring file variants against one chromosome's shards.

        No `min_overlap` here: the match rate is genome-wide, so only
        MATCH_COMBINE can decide it (docs/tickets/0003).

        Returns the chromosome that was matched.
        """
        logger.info("Matching scoring file variants in target genomes")

        with self._connect() as conn:
            # read variant information from zarr arrays into DB
            self._load_variants_from_zarr(conn)

            # check the shards can be scored independently before doing the work
            chromosome = check_cache_shape(
                conn=conn, sampleset_group=self.sampleset_group
            )

            # match score variants to zarr variants
            update_match_table(
                conn=conn,
                match_ambiguous=match_ambiguous,
                match_multiallelic=match_multiallelic,
                sampleset=self.sampleset,
            )
            return chromosome

    @property
    def sampleset_group(self) -> zarr.Group:
        return cast("zarr.Group", self.genotypes_root[self.sampleset])

    @property
    def shard(self) -> str:
        """The single leaf group this SCORE task holds.

        weight_mat_row_nr indexes this task's dosage array, so staging two shards
        at once would silently misalign it (docs/tickets/0007).
        """
        shards = list(get_shard_attrs(self.sampleset_group))
        if len(shards) != 1:
            raise ValueError(
                f"SCORE runs per shard but this cache holds {sorted(shards)}. "
                "Stage one shard per task."
            )
        return shards[0]

    def export_match_table(self, out_path: Pathish) -> None:
        """Write this chromosome's match table for MATCH_COMBINE to pool"""
        logger.info(f"Exporting match table to {str(out_path)}")
        with self._connect() as conn:
            conn.execute(
                "COPY allele_match_table TO ? (FORMAT parquet, COMPRESSION zstd)",
                [str(out_path)],
            )

    def export_samples(self, out_path: Pathish, chromosome: str) -> None:
        """Write the chromosome's sample IDs so MATCH_COMBINE can check the union"""
        logger.info(f"Exporting samples to {str(out_path)}")
        samples = sorted(
            sample
            for attrs in get_shard_attrs(self.sampleset_group).values()
            for sample in attrs["samples"]
        )
        with self._connect() as conn:
            conn.execute(
                """
                COPY (
                    SELECT $sampleset AS sampleset,
                           $chromosome AS chromosome,
                           UNNEST($samples) AS sample_id
                ) TO $path (FORMAT parquet, COMPRESSION zstd)
                """,
                {
                    "sampleset": self.sampleset,
                    "chromosome": chromosome,
                    "samples": samples,
                    "path": str(out_path),
                },
            )

    def import_match_tables(self, paths: Sequence[Pathish]) -> None:
        """Load partial match tables written by MATCH (or the gated table written
        by MATCH_COMBINE, which drops the target allele columns)"""
        logger.info(f"Importing {len(paths)} match tables")
        with self._connect() as conn:
            # score_log_table declares enum-typed columns, and update_match_table
            # (the only other caller) never runs in a combine task
            _setup_enums(conn)
            create_allele_match_table(conn)
            conn.execute(
                """
                INSERT INTO allele_match_table BY NAME
                SELECT * FROM read_parquet($paths)
                """,
                {"paths": [str(p) for p in paths]},
            )

    def combine_matches(self, min_overlap: float, out_path: Pathish) -> None:
        """Pool every chromosome's matches, apply min_overlap, and write the
        genome-wide match table restricted to accessions that passed"""
        if not (0 < min_overlap <= 1):
            raise ValueError(
                f"min_overlap must be positive and less than 1, got {min_overlap}"
            )

        with self._connect() as conn:
            add_match_log(
                conn=conn,
                min_overlap=min_overlap,
                sampleset=self.sampleset,
            )

        accessions = get_ok_accessions(db_path=self.db_path, sampleset=self.sampleset)

        logger.info(f"Exporting gated match table to {str(out_path)}")
        with self._connect() as conn:
            conn.execute(
                """
                COPY (
                    SELECT * EXCLUDE (target_ref, target_alts)
                    FROM allele_match_table
                    WHERE accession IN $accessions
                ) TO $path (FORMAT parquet, COMPRESSION zstd)
                """,
                {"accessions": accessions, "path": str(out_path)},
            )

    def export_partial_scores(self, out_path: Pathish) -> None:
        """Write this shard's contribution for SCORE_COMBINE to sum.

        Columns are listed rather than starred to pin the artefact schema, and
        `score_avg` is left out because it is GENERATED and regenerated on the
        combined table. `zarr_group` is what lets SCORE_COMBINE re-assert the
        uniqueness constraint that parquet cannot carry.
        """
        logger.info(f"Exporting partial scores to {str(out_path)}")
        with self._connect() as conn:
            conn.execute(
                """
                COPY (
                    SELECT sampleset, accession, n_matched, sample_id,
                           allele_count, score, dosage_sum,
                           $zarr_group AS zarr_group
                    FROM score_table
                ) TO $path (FORMAT parquet, COMPRESSION zstd)
                """,
                {"zarr_group": self.shard, "path": str(out_path)},
            )

    def combine_scores(self, paths: Sequence[Pathish]) -> None:
        """Sum every shard's partial scores.

        Valid without any awareness of how the cache was split, because every
        (sample, variant) pair appears in exactly one shard (docs/tickets/0008).
        """
        create_score_table(db_path=self.db_path)

        logger.info(f"Combining {len(paths)} partial score tables")
        parquet_paths = [str(p) for p in paths]

        with self._connect() as conn:
            # n_rows / n_shards ride along on the aggregate that has to happen
            # anyway, so checking costs a pass over the combined table rather than
            # a second pass over the partials
            conn.execute(
                """
                CREATE OR REPLACE TEMP TABLE combined AS
                SELECT sampleset,
                       accession,
                       sample_id,
                       -- already genome-wide (COUNT(DISTINCT row_nr) over the
                       -- pooled match table), so carried through, never summed
                       any_value(n_matched)       AS n_matched,
                       SUM(allele_count)          AS allele_count,
                       SUM(score)                 AS score,
                       SUM(dosage_sum)            AS dosage_sum,
                       COUNT(*)                   AS n_rows,
                       COUNT(DISTINCT zarr_group) AS n_shards
                FROM read_parquet($paths)
                GROUP BY sampleset, accession, sample_id
                """,
                {"paths": parquet_paths},
            )

            check_shards_contribute_once(conn=conn, paths=parquet_paths)

            # column order follows the score_table DDL; score_avg is generated
            conn.execute("""
                INSERT INTO score_table
                SELECT sampleset, accession, n_matched, sample_id,
                       allele_count, score, dosage_sum
                FROM combined
                """)

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

            # the match table was gated by MATCH_COMBINE, so every accession left
            # in it passed min_overlap and none of the rejected ones cost compute.
            # Scoped to this shard: an accession with no variants here contributes
            # nothing, and the samples it covers are scored by another shard
            # (every chromosome covers the same samples - check_sample_union)
            with duckdb.connect(str(self.db_path), read_only=True) as conn:
                accessions_to_calculate = [
                    row[0]
                    for row in conn.execute(
                        """
                        SELECT DISTINCT accession FROM allele_match_table
                        WHERE sampleset = ? AND zarr_group = ?
                          AND match_result.match_summary = 'matched'
                        ORDER BY accession
                        """,
                        [sampleset, self.shard],
                    ).fetchall()
                ]

            if not accessions_to_calculate:
                logger.info(
                    f"No matched variants in {self.shard!r}; exporting no contributions"
                )
                continue

            # get wide effect weights (a matrix) for each zarr genotype group
            grouped_dfs: dict[str, pl.DataFrame] = store_group_weight_arrays(
                db_path=self.db_path,
                sampleset=sampleset,
                pgs_group=pgs_group,
                accessions=accessions_to_calculate,
                zarr_group=self.shard,
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


def get_shard_attrs(sampleset_group: zarr.Group) -> dict[str, dict[str, Any]]:
    """Map every shard's path to its zarr attributes.

    Shards are the only nodes that carry sample IDs, so that's how they're found -
    without assuming how deep the hierarchy is.
    """
    return {
        node.path: dict(node.attrs)
        for _, node in sampleset_group.members(max_depth=None)
        if "samples" in node.attrs
    }


def check_cache_shape(
    conn: duckdb.DuckDBPyConnection, sampleset_group: zarr.Group
) -> str:
    """Check the shards under this chromosome can be scored independently.

    Reads zarr attributes and variant metadata only - nothing touches the 3D
    genotype arrays. Returns the chromosome (docs/tickets/0005).
    """
    shards = get_shard_attrs(sampleset_group)

    if not shards:
        raise ValueError(f"No shards found under {sampleset_group.path!r}")

    # leaves are {sampleset}/{chromosome}/{filename}, so a MATCH task's cache must
    # hold exactly one chromosome
    chromosomes = {path.split("/")[1] for path in shards}
    if len(chromosomes) != 1:
        raise ValueError(
            f"MATCH runs per chromosome but this cache holds {sorted(chromosomes)}. "
            "Stage one chromosome's shards per task."
        )

    check_absence_semantics(shards)
    check_samples_disjoint(shards)
    check_shards_are_rectangular(conn, shards)

    return chromosomes.pop()


def check_absence_semantics(shards: dict[str, dict[str, Any]]) -> None:
    """Absence semantics is a property of the callset, so every shard of it must
    declare the same thing, and only `unobserved` is implemented (ADR 0004)"""
    declared = {
        path: attrs.get("absence_semantics", "unobserved")
        for path, attrs in shards.items()
    }

    if len(set(declared.values())) > 1:
        raise ValueError(
            f"Shards disagree about absence semantics: {declared}. Shards of one "
            "sampleset are slices of one callset, so they cannot mean different "
            "things by an absent variant."
        )

    if (semantics := next(iter(declared.values()))) != "unobserved":
        raise NotImplementedError(
            f"absence_semantics={semantics!r} is not implemented yet. Absent "
            "variants are always treated as never assessed in this release."
        )


def check_samples_disjoint(shards: dict[str, dict[str, Any]]) -> None:
    """Shards under one chromosome must not share samples.

    SCORE_COMBINE sums partial scores grouped by sample ID, which is only valid
    while every (sample, variant) pair appears in exactly one shard. A shared
    sample double-counts its dosage: no error, no NaN, just an inflated score.
    """
    # ponytail: O(n^2) over shards-per-chromosome, which is a handful of sample
    # batches. If that ever stops being true, count sample occurrences in one pass.
    samples = {path: set(attrs["samples"]) for path, attrs in shards.items()}

    for (left, left_samples), (right, right_samples) in itertools.combinations(
        samples.items(), 2
    ):
        if shared := left_samples & right_samples:
            examples = ", ".join(sorted(shared)[:3])
            raise ValueError(
                f"Shards {left!r} and {right!r} share {len(shared)} samples "
                f"(e.g. {examples}). Shards under one chromosome must hold disjoint "
                "samples: the chromosome is the variant axis, so splitting variants "
                "within one is not supported. Split by chromosome instead."
            )


def check_shards_are_rectangular(
    conn: duckdb.DuckDBPyConnection, shards: dict[str, dict[str, Any]]
) -> None:
    """Shards under one chromosome must agree on which variants they contain.

    Defined over rows with real data: the cache back-fills queried-but-absent
    positions as sentinel rows with a null ref, so raw position sets are identical
    across shards by construction and would never look ragged.
    """
    if len(shards) < 2:
        return

    # targetvariants is keyed on (variant_id, zarr_group), so a variant can appear at
    # most once per shard and counting its rows counts the shards holding it. Anything
    # short of every shard is a variant one shard assessed and another didn't.
    ragged = conn.execute(
        """
        SELECT variant_id, list(zarr_group ORDER BY zarr_group) AS present_in
        FROM targetvariants
        WHERE ref IS NOT NULL AND zarr_group IN $shards
        GROUP BY variant_id
        HAVING COUNT(*) <> $n_shards
        ORDER BY variant_id
        LIMIT 3
        """,
        {"shards": list(shards), "n_shards": len(shards)},
    ).fetchall()

    if not ragged:
        return

    # only on the error path, so the happy path pays for one aggregate and no more
    counts = conn.execute(
        """
        SELECT zarr_group, COUNT(*)
        FROM targetvariants
        WHERE ref IS NOT NULL AND zarr_group IN $shards
        GROUP BY zarr_group
        ORDER BY zarr_group
        """,
        {"shards": list(shards)},
    ).fetchall()

    summary = ", ".join(f"{path} ({n} variants)" for path, n in counts)
    examples = ", ".join(
        f"{variant_id} (only in {', '.join(present_in)})"
        for variant_id, present_in in ragged
    )
    raise ValueError(
        "Shards under one chromosome contain different variants: "
        f"{summary}. For example: {examples}. Every shard of a chromosome must "
        "interrogate the same variants, because a variant missing from one shard and "
        "present in another cannot be told apart from a variant that was never "
        "assessed there. Re-export the target genomes without dropping sites that are "
        "monomorphic within a batch, or split the input by chromosome rather than by "
        "sample."
    )


def check_shards_contribute_once(
    conn: duckdb.DuckDBPyConnection, paths: Sequence[Pathish]
) -> None:
    """Each shard must reach SCORE_COMBINE exactly once.

    Reads the `combined` temp table built by `combine_scores`, where a sample's
    partial scores have already been counted: they must come from as many distinct
    shards as there are rows to sum.

    This is the constraint the `score_table` primary key cannot express. Within a
    shard the key already guarantees one row per (accession, sample), but a
    duplicate arrives as a separate file - a glob collecting the same output twice,
    or a retried task leaving its earlier output alongside the new one - and parquet
    carries no keys across that boundary. Summing it twice doubles that shard's
    score, allele count and dosage sum with no error and no NaN.
    """
    duplicated = conn.execute("""
        SELECT accession, sample_id, n_rows, n_shards
        FROM combined
        WHERE n_rows <> n_shards
        ORDER BY accession, sample_id
        LIMIT 3
        """).fetchall()

    if not duplicated:
        return

    # only on the error path, so the happy path pays nothing to name the shard
    accession, sample_id, _, _ = duplicated[0]
    shards = conn.execute(
        """
        SELECT zarr_group, COUNT(*) AS n
        FROM read_parquet($paths)
        WHERE accession = $accession AND sample_id = $sample_id
        GROUP BY zarr_group
        HAVING COUNT(*) > 1
        ORDER BY zarr_group
        """,
        {
            "paths": [str(p) for p in paths],
            "accession": accession,
            "sample_id": sample_id,
        },
    ).fetchall()

    listed = ", ".join(f"{group} ({n} times)" for group, n in shards)
    affected = ", ".join(f"{acc}/{sample}" for acc, sample, _, _ in duplicated)
    raise ValueError(
        f"Shard {listed} reached SCORE_COMBINE more than once. Partial scores are "
        "summed per sample, so scoring a shard twice doubles its score, allele "
        f"count and dosage sum. Affected (accession/sample): {affected}. Check that "
        "each shard contributes exactly one partial_scores.parquet - a stale output "
        "from a retried task, or a glob matching the same file twice, both do this."
    )


def check_sample_union(
    conn: duckdb.DuckDBPyConnection, paths: Sequence[Pathish]
) -> None:
    """Every chromosome in the cache must cover the same samples.

    Scoped to chromosomes present in the cache, which is automatic: a chromosome
    whose index query returned nothing is deleted by `cleanup_empty_zarr_group`,
    runs no MATCH task and writes no samples file (docs/tickets/0004).
    """
    mismatched = conn.execute(
        """
        WITH samples AS (SELECT * FROM read_parquet($paths)),
             union_of_all AS (SELECT DISTINCT sampleset, sample_id FROM samples),
             chromosomes AS (SELECT DISTINCT sampleset, chromosome FROM samples),
             expected AS (SELECT * FROM chromosomes JOIN union_of_all
                          USING (sampleset))
        SELECT DISTINCT chromosome
        FROM (SELECT * FROM expected EXCEPT SELECT * FROM samples)
        """,
        {"paths": [str(p) for p in paths]},
    ).fetchall()

    if mismatched:
        raise ValueError(
            f"Chromosomes {sorted(x[0] for x in mismatched)} are missing samples "
            "that other chromosomes have. Every chromosome of a sampleset must "
            "cover the same samples, or scores would be summed over a different "
            "cohort per chromosome."
        )


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
            # the path is relative to the store root, so it identifies the shard
            # whatever the depth of the hierarchy (e.g. UKB/1/batch1.vcf.gz)
            yield group.path, meta
        else:
            # recurse the group tree
            yield from get_variants_from_zarr(subgrp)
