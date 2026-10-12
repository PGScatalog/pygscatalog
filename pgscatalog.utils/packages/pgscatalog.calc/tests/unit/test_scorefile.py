import pathlib

import duckdb
import numpy as np
import polars as pl
import pytest
import zarr


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

    query = (
        duckdb.read_csv([pgs001229, pgs000001], dtype={"chr_name": str})
        .select("chr_name", "chr_position")
        .filter("chr_name is not null and chr_position is not null")
    )

    query = query.filter(f"chr_name IN {list(VALID_CHROMOSOMES - {None})}")

    assert scorefile_variants == frozenset(query.distinct().fetchall())

    # test optional chrom filter
    scorefile_variants_chr1 = frozenset(sfs.get_unique_positions(chrom="1"))
    query_chr1 = (
        duckdb.read_csv([pgs001229, pgs000001], dtype={"chr_name": str})
        .select("chr_name", "chr_position")
        .filter("chr_name is not null and chr_position is not null")
        .filter("chr_name = '1'")
    )

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
    meta.create_array(name="chr_name", data=np.array(cached_chr, dtype="U"))
    meta.create_array(name="chr_pos", data=np.array(cached_pos, dtype="int64"))

    # scorefile anti join with cached positions, should exclude cached positions
    uncached = one_sf.get_unique_positions(zarr_group=zarr_group)

    assert set(uncached).union(set(cached)) == set(all_variants)


def test_load_scoring_files(pgs000001, pgs001229):
    from pgscatalog.calc.lib.scorefile import load_scoring_files

    # an in-memory connection is enough, so no database on disk
    with duckdb.connect(":memory:") as conn:
        load_scoring_files(conn=conn, scorefile_paths=[pgs000001, pgs001229])

        # two unique accessions should be in the db table
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


@pytest.fixture
def xy_scorefile(tmp_path) -> str:
    """One XY variant and one X variant at the same position, plus a chr-prefixed
    autosome"""
    path = tmp_path / "PGS999997.txt"
    path.write_text(
        "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\t"
        "effect_type\tis_duplicated\taccession\trow_nr\n"
        "XY\t2781479\tA\tG\t0.1\tadditive\tFalse\tPGS999997\t0\n"
        "X\t2781479\tA\tG\t0.2\tadditive\tFalse\tPGS999997\t1\n"
        "chr7\t1000\tA\tG\t0.3\tadditive\tFalse\tPGS999997\t2\n"
    )
    return str(path)


def test_scoring_file_chromosomes_are_canonicalised(xy_scorefile):
    """The scoring file loader is a boundary, so MATCH only ever sees canonical
    names (ADR 0008)"""
    from pgscatalog.calc.lib.scorefile import load_scoring_files

    with duckdb.connect(":memory:") as conn:
        load_scoring_files(conn=conn, scorefile_paths=[xy_scorefile])
        # documented as re-runnable, and the UDF registration mustn't break that
        load_scoring_files(conn=conn, scorefile_paths=[xy_scorefile])
        rows = conn.sql(
            """
            SELECT chr_name, chr_position, is_duplicated, row_nr
            FROM score_variant_table ORDER BY row_nr
            """
        ).fetchall()

    # recorded behaviour, not a policy (ADR 0008 leaves this open): XY and X at one
    # position become two rows at one canonical position. Nothing collapses them and
    # nothing raises, and is_duplicated stays whatever the scoring file said, because
    # pgscatalog-format computed it before XY became X
    assert rows == [
        ("X", 2781479, False, 0),
        ("X", 2781479, False, 1),
        ("7", 1000, False, 2),
    ]


def test_unique_positions_collapse_xy_into_x(xy_scorefile):
    """Position queries are DISTINCT, so the target genome is only asked once"""
    from pgscatalog.calc import Scorefiles

    positions = Scorefiles(xy_scorefile).get_unique_positions(chrom="X")

    assert positions == [("X", 2781479)]


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
        {"chr_name": "2", "chr_pos": 30000},
    ]


def test_position_df_without_meta(tmp_path):
    from pgscatalog.calc.lib.scorefile import get_position_df

    zarr_group = zarr.group()
    df = get_position_df(zarr_group)

    assert isinstance(df, pl.DataFrame)
    assert df.columns == ["chr_name", "chr_pos"]
    assert df.height == 0
    assert df.to_dicts() == []
