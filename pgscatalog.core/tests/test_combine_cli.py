import collections
import csv
import gzip
import itertools
import json
from unittest.mock import patch
import pytest

from pgscatalog.core.cli.combine_cli import run
from pgscatalog.core import ScoringFile


@pytest.fixture(scope="package")
def custom_scorefiles(request):
    custom1 = request.path.parent / "data" / "custom.txt"
    custom2 = request.path.parent / "data" / "custom_no_oa.txt"
    return custom1, custom2


@pytest.fixture(scope="package", params=("GRCh37", "GRCh38"))
def harmonised_scorefiles(request):
    pgs000001 = request.path.parent / "data" / f"PGS000001_hmPOS_{request.param}.txt.gz"
    pgs000002 = request.path.parent / "data" / f"PGS000002_hmPOS_{request.param}.txt.gz"
    return (request.param, pgs000001), (request.param, pgs000002)


@pytest.fixture(scope="package")
def non_additive_scorefile_grch38(request):
    return request.path.parent / "data" / "PGS004255_hmPOS_GRCh38.txt.gz"


@pytest.fixture(scope="package")
def pgs000001_grch38(request):
    return request.path.parent / "data" / "PGS000001_hmPOS_GRCh38.txt.gz"


@pytest.fixture(scope="package", params=("GRCh37", "GRCh38"))
def scorefiles(request):
    pgs000001 = request.path.parent / "data" / "PGS000001.txt.gz"
    pgs000002 = request.path.parent / "data" / "PGS000002.txt.gz"
    return (request.param, pgs000001), (request.param, pgs000002)


@pytest.fixture(scope="package")
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


@pytest.fixture(scope="package")
def n_variants():
    return collections.Counter({"PGS000001": 77, "PGS000002": 77})


@pytest.fixture(scope="package")
def chain_dir(request):
    return request.path.parent / "data" / "chain"


@pytest.fixture(scope="package")
def lift_scorefiles(request):
    """These scoring files had chr_position column added for liftover test,
    genome_build modified in the header, and harmonisation columns deleted.
    The tuple is (target_build, score_path)"""
    return (
        (
            "GRCh38",
            request.path.parent / "data" / "lift" / "PGS000001_hmPOS_GRCh37.txt",
        ),
        (
            "GRCh37",
            request.path.parent / "data" / "lift" / "PGS000001_hmPOS_GRCh38.txt",
        ),
    )


def test_combine_nonadditive(tmp_path, non_additive_scorefile_grch38):
    """Test normalising a single non-additive scoring file fails."""
    out_path = tmp_path / "combined.txt.gz"
    path = [str(non_additive_scorefile_grch38)]

    args = [("pgscatalog-combine", "-s"), path, ("-o", str(out_path), "-t", "GRCh38")]
    flargs = list(itertools.chain(*args))

    with pytest.raises(ValueError):
        with patch("sys.argv", flargs):
            run()


def test_combine_skip(
    tmp_path, non_additive_scorefile_grch38, pgs000001_grch38, n_variants
):
    """Test that combining skips non-additive files but otherwise completes successfully."""
    out_path = tmp_path / "combined.txt.gz"
    paths = [str(non_additive_scorefile_grch38), str(pgs000001_grch38)]

    args = [("pgscatalog-combine", "-s"), paths, ("-o", str(out_path), "-t", "GRCh38")]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run()

    with gzip.open(out_path, mode="rt") as f:
        csv_reader = csv.DictReader(f, delimiter="\t")
        results = list(csv_reader)

    # split to remove harmonisation suffix hmPOS_GRCh3X
    pgs = collections.Counter([x["accession"].split("_")[0] for x in results])
    assert pgs["PGS000001"] == n_variants["PGS000001"]

    # the log should contain records of two scoring files though:
    with open(tmp_path / "log_combined.json") as f:
        log = json.load(f)
        assert len(log) == 2, "Missing scorefile from log"
        assert sum(x["compatible_effect_type"] is True for x in log) == 1
        assert sum(x["compatible_effect_type"] is False for x in log) == 1


def test_combine_score(tmp_path, scorefiles, expected_fields, n_variants):
    """Test combining unharmonised files.
    Genome build is missing from these files, so it should fail."""
    out_path = tmp_path / "combined.txt.gz"
    paths = [str(x) for _, x in scorefiles]
    build = "GRCh38"

    args = [("pgscatalog-combine", "-s"), paths, ("-o", str(out_path), "-t", build)]
    flargs = list(itertools.chain(*args))

    with pytest.raises(ValueError) as excinfo:
        with patch("sys.argv", flargs):
            run()
        assert "Can't combine files with missing build" in str(excinfo.value)


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

    # split to remove harmonisation suffix hmPOS_GRCh3X
    pgs = collections.Counter([x["accession"].split("_")[0] for x in results])
    assert pgs == n_variants
    assert all([expected_fields == tuple(variant.keys()) for variant in results])

    assert (tmp_path / "log_combined.json").exists(), "Missing log output"


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

    with pytest.raises(ValueError) as excinfo:
        with patch("sys.argv", flargs):
            run()
        assert "without --liftover" in str(excinfo.value)


def test_combine_custom(tmp_path, custom_scorefiles):
    """Test combining custom scoring files (not from PGS Catalog), including a scoring file that's missing the other allele column"""
    out_path = tmp_path / "combined.txt"
    target_build = "GRCh37"

    args = [
        ("pgscatalog-combine", "-s"),
        (str(x) for x in custom_scorefiles),
        ("-o", str(out_path), "-t", target_build),
    ]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run()

    with open(out_path, mode="rt") as f:
        csv_reader = csv.DictReader(f, delimiter="\t")
        results = list(csv_reader)

    assert results == [
        {
            "chr_name": "1",
            "chr_position": "10397",
            "effect_allele": "C",
            "other_allele": "A",
            "effect_weight": "-2.76E-02",
            "effect_type": "additive",
            "is_duplicated": "False",
            "accession": "test",
            "row_nr": "0",
        },
        {
            "chr_name": "1",
            "chr_position": "10397",
            "effect_allele": "C",
            "other_allele": "",
            "effect_weight": "-2.76E-02",
            "effect_type": "additive",
            "is_duplicated": "False",
            "accession": "test2",
            "row_nr": "0",
        },
    ]


def test_liftover(tmp_path, chain_dir, lift_scorefiles):
    """Test lifting over a scoring file from GRCh37 to GRCh38.
    Compare lifted coordinates with known good results"""
    out_path = tmp_path / "combined.txt"
    target_build, path = lift_scorefiles[0]
    _, grch38_path = lift_scorefiles[1]
    args = [
        ("pgscatalog-combine", "-s", str(path)),
        ("-o", str(out_path), "-t", target_build),
        ("--liftover", "--chain_dir", str(chain_dir)),
    ]

    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run()

    with open(out_path) as f:
        lifted_reader = csv.DictReader(f, delimiter="\t")
        known_good = [x.chr_position for x in ScoringFile(grch38_path).variants]
        lifted_pos = [row["chr_position"] for row in lifted_reader]

    assert all([int(x) == int(y) for x, y in zip(known_good, lifted_pos, strict=True)])
