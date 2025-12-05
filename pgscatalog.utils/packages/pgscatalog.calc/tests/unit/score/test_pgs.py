import random

import duckdb
import numpy as np
import polars as pl
import pytest
import zarr

from pgscatalog.calc.lib.score._pgs import (
    calculate_scores, 
    insert_scores,
    create_score_table, 
    write_scores,
    calculate_score_statistics,
    calculate_effect_allele_dosage,
    calculate_nonmissing_allele_count,
    )


@pytest.fixture
def tiny_dosage_matrix():
    """
    4 variants x 3 samples ->dosage matrix.
    3 variants belongs to PGS000001, while 1 in the PGS999999

    Variants: v0, v1, v2, v3
    Samples:  sample1, sample2, sample3

    v0: [0, 1, 2] -> PGS000001
    v1: [1, 0, 1] -> PGS000001
    v2: [2, 1, 0] -> PGS000001
    v3: [2, 2, 2] -> PGS999999

    Suppose everyone has ploidy of 2
    """
    return np.array(
        [
            [0.0, 1.0, 2.0],  # v0
            [1.0, 0.0, 1.0],  # v1
            [2.0, 1.0, 0.0],  # v2
            [2.0, 2.0, 2.0],  # v3
        ],
        dtype=np.float64,
        )


@pytest.fixture
def tiny_is_missing_matrix():
    """
    Boolean mask of missing calls (True = missing) in the  tiny_dosage_matrix
    """
    return np.array(
        [
            [False, False, True],   # v0
            [False, True,  False],  # v1
            [True, False, False],  # v2
            [True,  False, False],  # v3
        ],
        dtype=bool,
        )


@pytest.fixture
def tiny_dosage_zarr(tiny_dosage_matrix):
    store = zarr.storage.MemoryStore()
    return zarr.create_array(store=store, data=tiny_dosage_matrix)


@pytest.fixture
def tiny_is_missing_zarr(tiny_is_missing_matrix):
    store = zarr.storage.MemoryStore()
    return zarr.create_array(store=store, data=tiny_is_missing_matrix)


@pytest.fixture
def tiny_sample_ids() -> list[str]:
    return ["sample1", "sample2", "sample3"]

@pytest.fixture
def tiny_effect_weights_zarr():
    """
    Simple effect weights for 4 variants (variants, 1) : all 1.0.
    """
    weights = np.ones((4, 1), dtype=np.float64)
    store = zarr.storage.MemoryStore()
    return zarr.create_array(store=store, data=weights)

@pytest.fixture
def tiny_effect_weights_zarr():
    """
    Simple effect weights for 4 variants (variants, 1) : all 1.0.
    """
    weights = np.ones((4, 1), dtype=np.float64)
    store = zarr.storage.MemoryStore()
    return zarr.create_array(store=store, data=weights)

@pytest.fixture
def zero_dosage():
    return np.zeros(shape=(4, 3)) # 4 variants, 3 samples


@pytest.fixture
def one_dosage():
    return np.ones(shape=(4, 3)) # 4 variants, 3 samples

@pytest.fixture
def constant_weights():
    return np.array([1, 1, 1, 1])[:, np.newaxis]

@pytest.fixture
def random_score_table_rows():
    d = {
        "sampleset": [],
        "accession": [],
        "n_matched": [],
        "sample_id": [],
        "allele_count": [],
        "score": [],
        "dosage_sum": [],
    }
    for i in range(5):
        d["sampleset"].append("test")
        d["accession"].append("PGStest")
        d["sample_id"].append(f"sample{i}")
        d["allele_count"].append(random.randint(1, 10))
        d["n_matched"].append(random.randint(1, 10))
        d["dosage_sum"].append(round(random.uniform(0.0, 2.0), 3))
        d["score"].append(round(random.uniform(1.0, 100.0), 2))

    return pl.DataFrame(d)

#----------- Test the function -------------------

def test_calculate_score_zero(zero_dosage, constant_weights, tmp_path):
    dosage_store = zarr.storage.MemoryStore()
    dosage_zarr = zarr.create_array(store=dosage_store, data=zero_dosage)

    weights_store = zarr.storage.MemoryStore()
    weights_zarr = zarr.create_array(store=weights_store, data=constant_weights)

    score = calculate_scores(dosage_zarr, weights_zarr).compute()
    # no effect alleles present in any samples
    assert np.array_equal(score, np.zeros(shape=(zero_dosage.shape[1]))[:, np.newaxis])


def test_calculate_score(one_dosage, constant_weights, tmp_path):
    dosage_store = zarr.storage.MemoryStore()
    dosage_zarr = zarr.create_array(store=dosage_store, data=one_dosage)

    weights_store = zarr.storage.MemoryStore()
    weights_zarr = zarr.create_array(store=weights_store, data=constant_weights)

    score = calculate_scores(dosage_zarr, weights_zarr).compute()
    # one copy of effect allele present in all samples
    expected_score = np.array([len(constant_weights) * 1] * one_dosage.shape[1])[
        :, np.newaxis
    ]
    assert np.array_equal(score, expected_score)

def test_calculate_score_shape_mismatch(one_dosage):

    dosage_store = zarr.storage.MemoryStore()
    dosage_zarr = zarr.create_array(store=dosage_store, data=one_dosage)

    # mismatched weights
    bad_weights = np.ones((one_dosage.shape[0] - 1, 1), dtype=np.float64)
    weights_store = zarr.storage.MemoryStore()
    weights_zarr = zarr.create_array(store=weights_store, data=bad_weights)

    with pytest.raises(ValueError, match="dosage and weights must have same shape"):
        _ = calculate_scores(dosage_zarr, weights_zarr).compute()

def test_calculate_scores_tiny(tiny_dosage_zarr, tiny_effect_weights_zarr):
    """
    v0: [0, 1, 2]
    v1: [1, 0, 1]
    v2: [2, 1, 0]
    v3: [2, 2, 2]
    sum: [5, 4, 5]
    """
    scores = calculate_scores(tiny_dosage_zarr, tiny_effect_weights_zarr).compute()

    # expected per-sample sums as a column vector
    expected = np.array([[5.0], [4.0], [5.0]], dtype=np.float64)

    assert scores.shape == expected.shape
    assert scores.dtype == np.float64
    assert np.array_equal(scores, expected)

def test_create_score_table(tmp_path):
    db_path = tmp_path / "test.db"
    create_score_table(db_path)

    with duckdb.connect(str(db_path)) as conn:
        cols = conn.sql(
            "SELECT column_name, column_type FROM (DESCRIBE TABLE score_table);"
        ).fetchall()

    assert cols == [
        ("sampleset", "VARCHAR"),
        ("accession", "VARCHAR"),
        ("n_matched", "UINTEGER"),
        ("sample_id", "VARCHAR"),
        ("allele_count", "UINTEGER"),
        ("score", "DOUBLE"),
        ("dosage_sum", "DOUBLE"),
        ("score_avg", "DOUBLE"),
    ]

def test_tiny_calculate_effect_allele_dosage(tiny_dosage_zarr):
    """
    PGS000001
    v0: [0, 1, 2]
    v1: [1, 0, 1]
    v2: [2, 1, 0]
    ---
    sum:[3, 2, 3]

    PGS999999
    v3: [2, 2, 2]
    """
    accession_idx = np.array([0, 1, 2], dtype=np.int64)  # use v0 and v1

    # returns a per-sample sum of dosages for all variants in this accession
    expected = np.array([3.0, 2.0, 3.0], dtype=np.float64)
    
    result = calculate_effect_allele_dosage(
        dosage_array=tiny_dosage_zarr,
        accession_idx=accession_idx,
        )

    assert result.shape == expected.shape
    assert result.dtype == np.float64
    assert np.array_equal(result, expected)

def test_calculate_nonmissing_allele_count_tiny(tiny_dosage_zarr, tiny_is_missing_zarr):
    """
    missing:
    v0: [0, 1, np.na]
    v1: [1, np.na, 1]
    v2: [np.na, 1, 0]
    
    per_smaple_allele count:
    [2,2,2]

    v3: [np.na, 2, 2]
    """
    accession_idx = np.array([0, 1, 2], dtype=np.int64)

    # missing variants will have been imputed, so replace these positions with np.nan
    # return a per-sample sum of non-missing allele count in this accession

    expected = np.array([2, 2, 2], dtype=np.int64)

    result = calculate_nonmissing_allele_count(
        is_missing_array=tiny_is_missing_zarr,
        dosage_array=tiny_dosage_zarr,
        accession_idx=accession_idx,
        )

    assert result.shape == expected.shape
    assert np.array_equal(result, expected)

def test_calculate_score_statistics_tiny(
    tmp_path,
    tiny_dosage_zarr,
    tiny_is_missing_zarr,
    tiny_sample_ids,
):
    db_path = tmp_path / "tiny_scores.db"

    with duckdb.connect(str(db_path)) as conn:
        # 1) Raw table with plain match_summary
        conn.sql(
            """
            CREATE TABLE allele_match_raw (
                sampleset TEXT,
                accession TEXT,
                match_summary TEXT
            );
            """
        )
        conn.sql(
            """
            INSERT INTO allele_match_raw VALUES
                ('test', 'PGS000001', 'matched'),
                ('test', 'PGS999999', 'matched');
            """
        )

        # 2) View with struct field 'match_result.match_summary'
        conn.sql(
            """
            CREATE VIEW allele_match_table AS
            SELECT
                sampleset,
                accession,
                struct_pack(match_summary := match_summary) AS match_result
            FROM allele_match_raw;
            """
        )

        # 3) wide_score_variants defining membership of each variant
        conn.sql(
            """
            CREATE TABLE wide_score_variants (
                weight_mat_row_nr INTEGER,
                "PGS000001" DOUBLE,
                "PGS999999" DOUBLE
            );
            """
        )
        conn.sql(
            """
            INSERT INTO wide_score_variants VALUES
                (0, 1.0, 0.0),  -- v0 in PGS000001
                (1, 1.0, 0.0),  -- v1 in PGS000001
                (2, 1.0, 0.0),  -- v2 in PGS000001
                (3, 0.0, 1.0);  -- v3 in PGS999999
            """
        )

    accessions = ["PGS000001", "PGS999999"]

    # Score matrix shape: (n_samples, n_accessions) = (3, 2)
    # Columns in order ['PGS000001', 'PGS999999']
    scores = np.array(
        [
            [0.1, 0.2],  # sample1
            [0.3, 0.4],  # sample2
            [0.5, 0.6],  # sample3
        ],
        dtype=np.float64,
    )

    # Run the function under test
    result = calculate_score_statistics(
        db_path=db_path,
        accessions=accessions,
        dosage_array=tiny_dosage_zarr,
        is_missing_array=tiny_is_missing_zarr,
        sampleset="test",
        score=scores,
        sample_ids=tiny_sample_ids,
    )

    # --- Expected values ---

    # From earlier:
    # PGS000001 uses v0, v1, v2
    #   dosage_sum:   [3, 2, 3]
    #   allele_count: [2, 2, 2]
    #
    # PGS999999 uses v3
    #   dosage_sum:   [2, 2, 2]
    #   allele_count: [0, 1, 1]

    expected = pl.DataFrame(
        {
            "accession": ["PGS000001"] * 3 + ["PGS999999"] * 3,
            "dosage_sum": [3.0, 2.0, 3.0] + [2.0, 2.0, 2.0],
            "allele_count": [2, 2, 2] + [0, 1, 1],
            "n_matched": [1] * 6,  # one row per accession in allele_match_raw
            "sample_id": tiny_sample_ids + tiny_sample_ids,
            "score": [0.1, 0.3, 0.5] + [0.2, 0.4, 0.6],
            "sampleset": ["test"] * 6,
        }
    )

    result_sorted = result.sort(["accession", "sample_id"])
    expected_sorted = expected.sort(["accession", "sample_id"])

    assert result_sorted.columns == [
        "accession",
        "dosage_sum",
        "allele_count",
        "n_matched",
        "sample_id",
        "score",
        "sampleset",
    ]
    assert result_sorted.columns == expected_sorted.columns

    for col in result_sorted.columns:
        if result_sorted[col].dtype in (pl.Float32, pl.Float64):
            np.testing.assert_allclose(
                result_sorted[col].to_numpy(),
                expected_sorted[col].to_numpy(),
            )
        else:
            assert result_sorted[col].to_list() == expected_sorted[col].to_list()

def test_insert_scores_tiny(tmp_path, random_score_table_rows):
    db_path = tmp_path / "tiny_scores.db"
    create_score_table(db_path)

    # Insert using the function under test
    insert_scores(db_path=db_path, _score_df=random_score_table_rows)

    with duckdb.connect(str(db_path)) as conn:
        result = conn.sql(
            """
            SELECT sampleset, accession, n_matched, sample_id,
                   allele_count, score, dosage_sum
            FROM score_table
            ORDER BY accession, sample_id;
            """
        ).pl()

    # We expect exactly the same rows (but maybe different order; we sorted both)
    expected = random_score_table_rows.sort(["accession", "sample_id"])
    result_sorted = result.sort(["accession", "sample_id"])

    assert result_sorted.shape == expected.shape
    for col in expected.columns:
        if expected[col].dtype in (pl.Float64, pl.Float32):
            np.testing.assert_allclose(
                result_sorted[col].to_numpy(),
                expected[col].to_numpy(),
            )
        else:
            assert result_sorted[col].to_list() == expected[col].to_list()

def test_write_scores(tmp_path, tmp_path_factory, random_score_table_rows):
    db_path = tmp_path / "test.db"
    out_path = tmp_path_factory.mktemp("out")
    create_score_table(db_path)

    with duckdb.connect(db_path) as conn:
        conn.sql(
            """
            INSERT INTO score_table
            SELECT sampleset, 
                   accession, 
                   n_matched,
                   sample_id,
                   allele_count,
                   score,
                   dosage_sum
            FROM random_score_table_rows
            """
        )

    write_scores(db_path=db_path, out_dir=out_path, max_memory_gb="1GB", threads=1)

    # test hive partitioning
    out_csv = out_path / "sampleset=test" / "accession=PGStest" / "data_0.csv.gz"
    assert out_csv.exists()

    # test output structure
    df = duckdb.read_csv(str(out_csv)).pl()
    assert df.shape[0] == len(random_score_table_rows)
    assert df.schema == {
        "sampleset": pl.String,
        "accession": pl.String,
        "sample_id": pl.String,
        "n_matched": pl.Int64,
        "allele_count": pl.Int64,
        "dosage_sum": pl.Float64,
        "score": pl.Float64,
        "score_avg": pl.Float64,
    }