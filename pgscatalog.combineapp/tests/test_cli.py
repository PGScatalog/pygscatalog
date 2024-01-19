import collections
import csv
import gzip
import itertools
from unittest.mock import patch
import pytest

from pgscatalog.combineapp.cli import run


@pytest.fixture(scope="module", params=("GRCh37", "GRCh38"))
def harmonised_scorefiles(request):
    pgs000001 = (
        request.path.parent / "testdata" / f"PGS000001_hmPOS_{request.param}.txt.gz"
    )
    pgs000002 = (
        request.path.parent / "testdata" / f"PGS000002_hmPOS_{request.param}.txt.gz"
    )
    return (request.param, pgs000001), (request.param, pgs000002)


@pytest.fixture(scope="module", params=("GRCh37", "GRCh38"))
def scorefiles(request):
    pgs000001 = request.path.parent / "testdata" / "PGS000001.txt.gz"
    pgs000002 = request.path.parent / "testdata" / "PGS000002.txt.gz"
    return (request.param, pgs000001), (request.param, pgs000002)


@pytest.fixture(scope="module")
def expected_fields():
    return (
        "chr_name",
        "chr_position",
        "effect_allele",
        "other_allele",
        "effect_weight",
        "effect_type",
        "is_duplicated",
        "accession",
        "row_nr",
    )


@pytest.fixture(scope="module")
def n_variants():
    return collections.Counter({"PGS000001": 77, "PGS000002": 77})


def test_combine_score(tmp_path, scorefiles, expected_fields, n_variants):
    """Test combining unharmonised files.
    Genome build is missing from these files, so it should fail."""
    out_path = tmp_path / "combined.txt.gz"
    paths = [str(x) for _, x in scorefiles]
    build = [str(x) for x, _ in scorefiles][0]

    args = [("pgscatalog-combine", "-s"), paths, ("-o", str(out_path), "-t", build)]
    flargs = list(itertools.chain(*args))

    with pytest.raises(NotImplementedError):
        with patch("sys.argv", flargs):
            run()


def test_combine_score_harmonised(
    tmp_path, harmonised_scorefiles, expected_fields, n_variants
):
    """Test combining harmonised data: PGS000001 and PGS000002.
    The genome build always matches across both files, so combining should work."""
    out_path = tmp_path / "combined.txt.gz"
    paths = [str(x) for _, x in harmonised_scorefiles]
    build = [str(x) for x, _ in harmonised_scorefiles][0]

    args = [("pgscatalog-combine", "-s"), paths, ("-o", str(out_path), "-t", build)]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run()

    with gzip.open(out_path, mode="rt") as f:
        csv_reader = csv.DictReader(f, delimiter="\t")
        results = list(csv_reader)

    pgs = collections.Counter([x["accession"] for x in results])
    assert pgs == n_variants
    assert all([expected_fields == tuple(variant.keys()) for variant in results])


def test_combine_fail(tmp_path, harmonised_scorefiles):
    """Test combining files with the same build but the wrong target build.
    Local files are in GRCh38, but we say we want GRCh37 (and vice versa)"""
    out_path = tmp_path / "combined.txt.gz"
    paths = [str(x) for _, x in harmonised_scorefiles]
    build = [str(x) for x, _ in harmonised_scorefiles][0]

    # uh oh
    if build == "GRCh37":
        target_build = "GRCh38"
    else:
        target_build = "GRCh37"

    args = [
        ("pgscatalog-combine", "-s"),
        paths,
        ("-o", str(out_path), "-t", target_build),
    ]
    flargs = list(itertools.chain(*args))

    # TODO: implement liftover?
    with pytest.raises(NotImplementedError):
        with patch("sys.argv", flargs):
            run()
