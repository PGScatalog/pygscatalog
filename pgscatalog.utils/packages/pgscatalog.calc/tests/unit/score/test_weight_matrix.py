import duckdb
import numpy as np
import polars as pl
import pytest
import zarr

from pgscatalog.calc.lib.score._weight_matrix import (
    create_wide_weights_table,
    get_variant_metadata,
    get_weight_matrix,
    store_group_weight_arrays,
    store_results_in_zarr,
)
from pgscatalog.calc.lib.score.scoredb import ScoreDatabase

"""
The _weight_matrix functions take a connection to a database with at least two tables:
score_variant_table and allele_match_table. Most of these tests never need the
database on disk, so they take an in-memory connection (docs/tickets/0016).
"""


def _load_test_tables(conn, score_variant_table, allele_match_table) -> None:
    """Load the pre-defined score_variant_table / allele_match_table test data"""
    conn.execute(
        """
        CREATE TABLE score_variant_table AS
        SELECT *
        FROM read_json(?, format='nd');
        """,
        [str(score_variant_table)],
    )
    conn.execute(
        """
        CREATE TABLE allele_match_table AS
        SELECT *
        FROM read_json(?, format='nd');
        """,
        [str(allele_match_table)],
    )


@pytest.fixture
def db_with_scores(allele_match_table, score_variant_table):
    """An in-memory database with two PGS scorefiles loaded for testing."""
    with duckdb.connect(":memory:") as conn:
        _load_test_tables(conn, score_variant_table, allele_match_table)
        yield conn


@pytest.fixture
def scoredb_with_scores(
    tmp_path, allele_match_table, score_variant_table
) -> ScoreDatabase:
    """The same data, but on disk.

    store_group_weight_arrays pivots and then reads back what it pivoted, so it needs
    a writable connection and then a read-only one, and duckdb can't open an in-memory
    database read-only.
    """
    db = ScoreDatabase(tmp_path / "weights.duckdb", max_memory_gb=1, threads=1)
    with db.connect() as conn:
        _load_test_tables(conn, score_variant_table, allele_match_table)
    return db


@pytest.fixture
def empty_scoredb(tmp_path) -> ScoreDatabase:
    """
    A database with the correct schema but no rows, used to test
    error handling in store_group_weight_arrays.
    """
    db = ScoreDatabase(tmp_path / "empty.duckdb", max_memory_gb=1, threads=1)

    with db.connect() as conn:
        conn.execute(
            """
            CREATE TABLE score_variant_table (
                accession TEXT,
                effect_type TEXT,
                row_nr UBIGINT,
                effect_weight DOUBLE
            );
            """
        )
        conn.execute(
            """
            CREATE TABLE allele_match_table (
                sampleset TEXT,
                accession TEXT,
                row_nr UBIGINT,
                match_result STRUCT(effect_allele_idx UINT8),
                target_row_nr UBIGINT,
                zarr_group TEXT
            );
            """
        )

    return db


@pytest.fixture
def zarr_group(tmp_path) -> zarr.Group:
    store_path = tmp_path / "weights.zarr"
    return zarr.open_group(store_path, mode="w")


# ---
def test_create_wide_weights_table(db_with_scores):
    conn = db_with_scores
    create_wide_weights_table(
        conn,
        sampleset="test",
        accessions="test",
        zarr_group="test/ALL/test.vcf.gz",
    )

    tables = conn.execute(
        """
            SELECT table_name
            FROM information_schema.tables
            WHERE table_schema = 'main';
            """
    ).fetchall()

    table_names = {table[0] for table in tables}
    assert "wide_score_variants" in table_names

    row = conn.execute(
        """
        SELECT
           weight_mat_row_nr,
           zarr_group,
           is_recessive,
           is_dominant,
           target_row_nr,
           effect_allele_idx,
           test
        FROM wide_score_variants;
        """
    ).fetchone()

    (
        test_weight_mat_row_nr,
        test_zarr_group,
        test_is_recessive,
        test_is_dominant,
        test_target_row_nr,
        test_effect_allele_idx,
        test_value,
    ) = row

    assert "wide_score_variants" in table_names
    assert test_weight_mat_row_nr == 0
    assert test_zarr_group == "test/ALL/test.vcf.gz"
    assert test_is_recessive is False
    assert test_is_dominant is False
    assert test_target_row_nr == 0
    assert test_effect_allele_idx == 1
    assert test_value == pytest.approx(-0.008818077)


def test_get_variant_metadata(db_with_scores):
    conn = db_with_scores
    group_path = "test/ALL/test.vcf.gz"

    create_wide_weights_table(
        conn,
        sampleset="test",
        accessions="test",
        zarr_group=group_path,
    )

    df = get_variant_metadata(conn=conn, group_path=group_path)

    assert isinstance(df, pl.DataFrame)
    assert df.shape == (1, 5)  # (n_variants, n_columns)

    # Columns
    assert df.columns == [
        "weight_mat_row_nr",
        "effect_allele_idx",
        "target_row_nr",
        "is_recessive",
        "is_dominant",
    ]


def test_get_weight_matrix(db_with_scores):
    conn = db_with_scores
    group_path = "test/ALL/test.vcf.gz"

    create_wide_weights_table(
        conn,
        sampleset="test",
        accessions="test",
        zarr_group=group_path,
    )

    cols, mat = get_weight_matrix(conn=conn, group_path=group_path)

    # Accessions should match the pivoted column names
    assert cols == ["test"]

    # Shape: matrix
    assert mat.shape == (1, 1)  # (n_variants, n_accessions)
    expected = np.array([[-0.008818077]])
    np.testing.assert_allclose(mat, expected)


def test_store_results_in_zarr(db_with_scores, zarr_group):
    conn = db_with_scores
    group_path = "test/ALL/test.vcf.gz"

    create_wide_weights_table(
        conn,
        sampleset="test",
        accessions="test",
        zarr_group=group_path,
    )

    store_results_in_zarr(conn=conn, group_path=group_path, pgs_group=zarr_group)

    assert "weight_matrix" in zarr_group.array_keys()
    arr = zarr_group["weight_matrix"]
    assert isinstance(arr, zarr.Array)

    assert arr.shape == (1, 1)
    expected = np.array([[-0.008818077]])
    np.testing.assert_allclose(np.array(arr), expected)

    assert zarr_group.attrs["accessions"] == ["test"]


def test_store_group_weight_arrays_success(scoredb_with_scores, zarr_group):
    df = store_group_weight_arrays(
        db=scoredb_with_scores,
        sampleset="test",
        pgs_group=zarr_group,
        accessions=["test"],
        zarr_group="test/ALL/test.vcf.gz",
    )

    assert isinstance(df, pl.DataFrame)
    assert df.shape == (1, 5)  # (n_variants, n_columns)
    assert df["weight_mat_row_nr"].to_list() == [0]
    assert df["target_row_nr"].to_list() == [0]


SHARD_A = "test/1/a.vcf.gz"
SHARD_B = "test/8/b.vcf.gz"


@pytest.fixture
def db_two_shards():
    """A genome-wide match table covering two shards, as MATCH_COMBINE writes it.

    target_row_nr is shard-local, so both shards number their variants from 0. In
    shard A, two scoring variants match the same target row with different effect
    allele indices (the REFALT/ALTREF case).
    """
    # (accession, row_nr, target_row_nr, effect_allele_idx, zarr_group)
    matches = [
        ("PGS000001", 0, 0, 1, SHARD_A),
        ("PGS000001", 1, 1, 1, SHARD_A),
        ("PGS000002", 2, 1, 0, SHARD_A),
        ("PGS000002", 3, 2, 1, SHARD_A),
        ("PGS000001", 4, 0, 1, SHARD_B),
        ("PGS000002", 5, 1, 1, SHARD_B),
    ]

    with duckdb.connect(":memory:") as conn:
        conn.execute("""
            CREATE TABLE score_variant_table (
                accession TEXT, effect_type TEXT, row_nr UBIGINT, effect_weight DOUBLE
            );
        """)
        conn.executemany(
            "INSERT INTO score_variant_table VALUES (?, 'additive', ?, ?)",
            [(accession, row_nr, 0.1) for accession, row_nr, *_ in matches],
        )
        conn.execute("""
            CREATE TABLE allele_match_raw (
                sampleset TEXT,
                accession TEXT,
                row_nr UBIGINT,
                target_row_nr UBIGINT,
                effect_allele_idx UINT8,
                zarr_group TEXT
            );
        """)
        conn.executemany(
            "INSERT INTO allele_match_raw VALUES ('test', ?, ?, ?, ?, ?)", matches
        )
        # match_result is a struct in the real schema
        conn.execute("""
            CREATE VIEW allele_match_table AS
            SELECT sampleset, accession, row_nr,
                   struct_pack(effect_allele_idx := effect_allele_idx) AS match_result,
                   target_row_nr, zarr_group
            FROM allele_match_raw;
        """)

        yield conn


def test_weights_are_scoped_to_one_shard(db_two_shards):
    """SCORE stages the genome-wide match table but holds one shard's genotypes.

    Without the zarr_group predicate, weights get built for a group whose arrays
    aren't in this task's cache, and weight_mat_row_nr stops being a valid index
    into this task's dosage array.
    """
    create_wide_weights_table(
        db_two_shards,
        sampleset="test",
        accessions=["PGS000001", "PGS000002"],
        zarr_group=SHARD_A,
    )

    rows = db_two_shards.execute("""
        SELECT weight_mat_row_nr, zarr_group, target_row_nr, effect_allele_idx
        FROM wide_score_variants ORDER BY weight_mat_row_nr
    """).fetchall()

    # shard B's two matches are absent, and weight_mat_row_nr is the position of a
    # row in shard A's own ordering, counting from 0
    assert rows == [
        (0, SHARD_A, 0, 1),
        (1, SHARD_A, 1, 0),
        (2, SHARD_A, 1, 1),
        (3, SHARD_A, 2, 1),
    ]

    # get_variant_metadata must agree, because write_dosage_arrays indexes the
    # genotype array with it
    df = get_variant_metadata(conn=db_two_shards, group_path=SHARD_A)
    assert df["weight_mat_row_nr"].to_list() == [0, 1, 2, 3]
    assert df["effect_allele_idx"].to_list() == [1, 0, 1, 1]


def test_refalt_altref_matches_get_a_row_each(db_two_shards):
    """One position matched REFALT by one score and ALTREF by another needs two
    weight matrix rows, because each one takes a different effect allele dosage.

    effect_allele_idx is the tie-break that makes them separate rows and gives the
    pair a stable order. Collapsing them would sum two accessions' weights onto one
    dosage row; leaving them unordered would let get_weight_matrix disagree with the
    weight_mat_row_nr assigned here.
    """
    create_wide_weights_table(
        db_two_shards,
        sampleset="test",
        accessions=["PGS000001", "PGS000002"],
        zarr_group=SHARD_A,
    )

    df = get_variant_metadata(conn=db_two_shards, group_path=SHARD_A)

    # target row 1 was matched by PGS000001 as ALT and PGS000002 as REF
    tied = df.filter(pl.col("target_row_nr") == 1)
    assert tied["effect_allele_idx"].to_list() == [0, 1]

    # weight_mat_row_nr indexes the dosage array, so it must be the position of the
    # row in this ordering - including across the tie
    assert df["weight_mat_row_nr"].to_list() == list(range(df.height))

    cols, mat = get_weight_matrix(conn=db_two_shards, group_path=SHARD_A)
    assert cols == ["PGS000001", "PGS000002"]
    np.testing.assert_allclose(
        mat,
        [
            [0.1, 0.0],  # target 0, ALT: PGS000001 only
            [0.0, 0.1],  # target 1, REF: PGS000002 only
            [0.1, 0.0],  # target 1, ALT: PGS000001 only
            [0.0, 0.1],  # target 2, ALT: PGS000002 only
        ],
    )


def test_store_group_weight_arrays_raises_on_empty(empty_scoredb, zarr_group):
    """
    When no variants are available for the given sampleset/accessions,
    store_group_weight_arrays should raise a ValueError.
    """
    with pytest.raises(ValueError):
        store_group_weight_arrays(
            db=empty_scoredb,
            sampleset="unknown_sampleset",
            pgs_group=zarr_group,
            accessions=["PGS000001"],
            zarr_group="unknown/ALL/unknown.vcf.gz",
        )
