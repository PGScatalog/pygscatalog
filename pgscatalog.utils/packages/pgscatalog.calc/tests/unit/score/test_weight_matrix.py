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

"""
The input file of the _weight_matrix is a path of a db with at least two tables: score_variant_table and allele_match_table.
"""


@pytest.fixture
def db_with_scores(tmp_path, allele_match_table, score_variant_table) -> str:
    """
    use the load_scoring_files function from the scorefile module to
    create a DuckDB file with two PGS scorefiles loaded for testing.
    """
    db_path = tmp_path / "weights.duckdb"
    conn = duckdb.connect(str(db_path))

    # Load the pre-defined test data score_variant_table
    conn.execute(
        """
        CREATE TABLE score_variant_table AS
        SELECT * 
        FROM read_json(?, format='nd');
        """,
        [str(score_variant_table)],
    )

    # Load the pre-defined test data allele_match_table
    conn.execute(
        """
        CREATE TABLE allele_match_table AS
        SELECT * 
        FROM read_json(?, format='nd');
        """,
        [str(allele_match_table)],
    )

    conn.close()
    return str(db_path)


@pytest.fixture
def empty_db(tmp_path) -> str:
    """
    DuckDB file with the correct schema but no rows, used to test
    error handling in store_group_weight_arrays.
    """
    db_path = tmp_path / "empty.duckdb"
    conn = duckdb.connect(str(db_path))

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
            filename TEXT
        );
        """
    )

    conn.close()
    return str(db_path)


@pytest.fixture
def zarr_group(tmp_path) -> zarr.Group:
    store_path = tmp_path / "weights.zarr"
    return zarr.open_group(store_path, mode="w")


# ---
def test_create_wide_weights_table(db_with_scores):
    db_path = db_with_scores
    create_wide_weights_table(db_path, sampleset="test", accessions="test")

    conn = duckdb.connect(str(db_path))

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
    assert test_zarr_group == "test/test.vcf.gz"
    assert test_is_recessive is False
    assert test_is_dominant is False
    assert test_target_row_nr == 0
    assert test_effect_allele_idx == 1
    assert test_value == pytest.approx(-0.008818077)


def test_get_variant_metadata(db_with_scores):
    db_path = db_with_scores

    # create the zarr_group
    group_path = create_wide_weights_table(
        db_path, sampleset="test", accessions="test"
    )[0]

    df = get_variant_metadata(db_path=db_with_scores, group_path=group_path)

    assert isinstance(df, pl.DataFrame)
    assert df.shape == (1, 6)  # (n_variants, n_columns)

    # Columns
    assert df.columns == [
        "zarr_group",
        "weight_mat_row_nr",
        "effect_allele_idx",
        "target_row_nr",
        "is_recessive",
        "is_dominant",
    ]


def test_get_weight_matrix(db_with_scores):
    db_path = db_with_scores

    # create the zarr_group
    group_path = create_wide_weights_table(
        db_path, sampleset="test", accessions="test"
    )[0]

    cols, mat = get_weight_matrix(db_path=db_with_scores, group_path=group_path)

    # Accessions should match the pivoted column names
    assert cols == ["test"]

    # Shape: matrix
    assert mat.shape == (1, 1)  # (n_variants, n_accessions)
    expected = np.array([[-0.008818077]])
    np.testing.assert_allclose(mat, expected)


def test_store_results_in_zarr(db_with_scores, zarr_group):
    db_path = db_with_scores

    # create the zarr_group
    group_path = create_wide_weights_table(
        db_path, sampleset="test", accessions="test"
    )[0]

    store_results_in_zarr(
        db_path=db_with_scores, group_path=group_path, pgs_group=zarr_group
    )

    assert "weight_matrix" in zarr_group.array_keys()
    arr = zarr_group["weight_matrix"]
    assert isinstance(arr, zarr.Array)

    assert arr.shape == (1, 1)
    expected = np.array([[-0.008818077]])
    np.testing.assert_allclose(np.array(arr), expected)

    assert zarr_group.attrs["accessions"] == ["test"]


def test_store_group_weight_arrays_success(db_with_scores, zarr_group):
    result = store_group_weight_arrays(
        db_path=db_with_scores,
        sampleset="test",
        pgs_group=zarr_group,
        accessions=["test"],
    )

    assert set(result.keys()) == {
        "test/test.vcf.gz",
    }

    df = result["test/test.vcf.gz"]

    assert isinstance(df, pl.DataFrame)
    assert df.shape == (1, 5)  # (n_variants, n_columns)
    assert df["weight_mat_row_nr"].to_list() == [0]
    assert df["target_row_nr"].to_list() == [0]


def test_store_group_weight_arrays_raises_on_empty(empty_db, zarr_group):
    """
    When no variants are available for the given sampleset/accessions,
    store_group_weight_arrays should raise a ValueError.
    """
    with pytest.raises(ValueError):
        store_group_weight_arrays(
            db_path=empty_db,
            sampleset="unknown_sampleset",
            pgs_group=zarr_group,
            accessions=["PGS000001"],
        )
