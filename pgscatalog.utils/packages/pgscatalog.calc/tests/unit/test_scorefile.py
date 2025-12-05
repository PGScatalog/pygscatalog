import pathlib

import duckdb
import pytest
import zarr
import polars as pl
import numpy as np


@pytest.fixture
def pgs001229() -> str:
    return str(
        pathlib.Path(__file__).parent.parent
        / "data"
        / "structured"
        / "PGS001229_hmPOS_GRCh38.txt.gz"
    )


@pytest.fixture
def pgs000001() -> str:
    return str(
        pathlib.Path(__file__).parent.parent
        / "data"
        / "structured"
        / "PGS000001_hmPOS_GRCh38.txt.gz"
    )


def test_scorefiles(tmp_path, pgs000001, pgs001229):
    from pgscatalog.calc import VALID_CHROMOSOMES, Scorefiles
    # test init with one scoring file
    one_sf = Scorefiles(pgs000001)
    assert [x.name for x in one_sf.paths] == ["PGS000001_hmPOS_GRCh38.txt.gz"]

    # test init with multiple scoring files
    sfs = Scorefiles([pgs000001, pgs001229])
    assert {x.name for x in sfs.paths} == {
        "PGS000001_hmPOS_GRCh38.txt.gz",
        "PGS001229_hmPOS_GRCh38.txt.gz",
    }

    # test getting unique variants manually
    scorefile_variants = frozenset(sfs.get_unique_positions())

    query = (duckdb.read_csv([pgs001229, pgs000001], dtype={"chr_name": str})
        .select("chr_name", "chr_position")
        .filter("chr_name is not null and chr_position is not null"))

    query = query.filter(
                        f"chr_name IN {list(VALID_CHROMOSOMES - {None})}"
                    )

    assert scorefile_variants == frozenset(query.distinct().fetchall())

    # test optional chrom filter
    scorefile_variants_chr1 = frozenset(sfs.get_unique_positions(chrom="1"))
    query_chr1 = (duckdb.read_csv([pgs001229, pgs000001], dtype={"chr_name": str})
        .select("chr_name", "chr_position")
        .filter("chr_name is not null and chr_position is not null")
        .filter("chr_name = '1'"))

    assert scorefile_variants_chr1 == frozenset(query_chr1.distinct().fetchall())

    # test when zarr group is provided
    all_variants = (
        duckdb.read_csv([pgs000001], dtype={"chr_name": str})
        .select("chr_name", "chr_position")
        .filter("chr_name is not null and chr_position is not null")
        .filter(f"chr_name IN {list(VALID_CHROMOSOMES - {None})}")
        .distinct()
    ).fetchall()

    cached = all_variants[:5]  # pretend these are already cached
    cached_chr = [c[0] for c in cached]
    cached_pos = [c[1] for c in cached]

    store = zarr.storage.MemoryStore()
    zarr_group = zarr.group(store=store, zarr_format=2)
    meta = zarr_group.create_group("meta")
    meta.create_array(name="chr_name",data=np.array(cached_chr, dtype="U"))
    meta.create_array(name="chr_pos", data=np.array(cached_pos, dtype="int64"))

    # scorefile anti join with cached positions, should exclude cached positions
    uncached = one_sf.get_unique_positions(zarr_group=zarr_group)

    assert set(uncached).union(set(cached)) == set(all_variants)

def test_load_scoring_files(tmp_path, pgs000001, pgs001229):
    from pgscatalog.calc.lib.scorefile import load_scoring_files

    db_path = tmp_path / "test.duckdb"

    load_scoring_files(
        db_path=db_path,
        scorefile_paths=[pgs000001, pgs001229],
        threads=1,
        max_memory_gb="4GB",
    )

    # two unique accessions should be in the db table
    with duckdb.connect(db_path) as conn:
        res = frozenset(
            conn.sql(
                """
                SELECT accession, COUNT(accession)
                FROM score_variant_table
                GROUP BY accession
                """
            ).fetchall()
        )

    assert res == frozenset([("PGS001229", 51183), ("PGS000001", 77)])

def test_position_df_with_meta(tmp_path):
    from pgscatalog.calc.lib.scorefile import get_position_df
    
    store = zarr.storage.MemoryStore()
    zarr_group = zarr.group(store=store, zarr_format=2)
    meta = zarr_group.create_group("meta")
    meta.create_array("chr_name", data=np.array(["1", "1", "2"], dtype="U"))
    meta.create_array("chr_pos", data=np.array([10000, 20000, 30000], dtype="int64"))

    df = get_position_df(zarr_group)

    assert df.columns == ["chr_name", "chr_pos"]
    assert df.to_dicts() == [
    {"chr_name": "1", "chr_pos": 10000},
    {"chr_name": "1", "chr_pos": 20000},
    {"chr_name": "2", "chr_pos": 30000}
    ]

def test_position_df_without_meta(tmp_path):
    from pgscatalog.calc.lib.scorefile import get_position_df
    
    zarr_group = zarr.group()
    df = get_position_df(zarr_group)

    assert isinstance(df, pl.DataFrame)
    assert df.columns == ["chr_name", "chr_pos"]
    assert df.height == 0
    assert df.to_dicts() == []
