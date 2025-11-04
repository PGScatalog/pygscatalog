import pathlib

import duckdb
import pytest

from pgscatalog.calc import Scorefiles


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


def test_scorefiles(pgs000001, pgs001229):
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

    variants = frozenset(
        duckdb.read_csv([pgs001229, pgs000001], dtype={"chr_name": str})
        .select("chr_name", "chr_position")
        .filter("chr_name is not null")
        .distinct()
        .fetchall()
    )

    assert scorefile_variants == variants
