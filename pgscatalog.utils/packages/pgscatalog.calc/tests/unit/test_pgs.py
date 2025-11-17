import random

import duckdb
import numpy as np
import polars as pl
import pytest
import zarr

from pgscatalog.calc.lib.score._pgs import calculate_scores, create_score_table, write_scores

RANDOM_SEED = 42


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


@pytest.fixture
def zero_dosage():
    return np.zeros(shape=(5, 1000))


@pytest.fixture
def one_dosage():
    return np.ones(shape=(5, 1000))


@pytest.fixture
def constant_weights():
    return np.array([1, 1, 1, 1, 1])[:, np.newaxis]


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
