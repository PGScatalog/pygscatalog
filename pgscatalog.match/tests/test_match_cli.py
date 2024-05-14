import itertools
import os
from unittest.mock import patch
import pytest

from pgscatalog.core import ZeroMatchesError, MatchRateError
from pgscatalog.match.cli.match_cli import run_match


@pytest.fixture(scope="module")
def good_scorefile(request):
    return request.path.parent / "data" / "good_match_scorefile.txt"


@pytest.fixture(scope="module")
def bad_scorefile(request):
    return request.path.parent / "data" / "bad_match_scorefile.txt.gz"


@pytest.fixture(scope="module")
def good_variants(request):
    return request.path.parent / "data" / "good_match.pvar"


def test_match(tmp_path_factory, good_scorefile, good_variants):
    """Test matching runs without errors with good data"""
    outdir = tmp_path_factory.mktemp("outdir")

    args = [
        (
            "pgscatalog-match",
            "-d",
            "test",
            "-s",
            str(good_scorefile),
            "-t",
            str(good_variants),
            "--outdir",
            str(outdir),
            "--min_overlap",
            "0.75",
        )
    ]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run_match()

    assert (outdir / "test_log.csv.gz").exists()
    assert (outdir / "test_summary.csv").exists()
    assert (outdir / "test_ALL_recessive_0.scorefile.gz").exists()
    assert (outdir / "test_ALL_dominant_0.scorefile.gz").exists()
    assert (outdir / "test_ALL_additive_0.scorefile.gz").exists()


def test_only_match(tmp_path_factory, good_scorefile, good_variants):
    """Test just matching (for big data)"""
    outdir = tmp_path_factory.mktemp("outdir")

    args = [
        (
            "pgscatalog-match",
            "-d",
            "test",
            "-s",
            str(good_scorefile),
            "-t",
            str(good_variants),
            "--outdir",
            str(outdir),
            "--only_match",
        )
    ]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run_match()

    # arrow IPC files have been written
    assert "0.ipc.zst" in os.listdir(outdir)


def test_strict_match(tmp_path_factory, good_scorefile, good_variants):
    """Test matching well with extremely strict overlap to trigger a MatchRateError"""
    outdir = tmp_path_factory.mktemp("outdir")

    args = [
        (
            "pgscatalog-match",
            "-d",
            "test",
            "-s",
            str(good_scorefile),
            "-t",
            str(good_variants),
            "--outdir",
            str(outdir),
            "--min_overlap",
            "1.0",
        )
    ]
    flargs = list(itertools.chain(*args))

    with pytest.raises(MatchRateError):
        with patch("sys.argv", flargs):
            run_match()

    assert (outdir / "test_log.csv.gz").exists()
    assert (outdir / "test_summary.csv").exists()
    assert not (outdir / "test_ALL_recessive_0.scorefile.gz").exists()
    assert not (outdir / "test_ALL_dominant_0.scorefile.gz").exists()
    assert not (outdir / "test_ALL_additive_0.scorefile.gz").exists()


def test_match_fail(tmp_path_factory, bad_scorefile, good_variants):
    """Test trying to match when no variants intersect"""
    outdir = tmp_path_factory.mktemp("outdir")

    args = [
        (
            "pgscatalog-match",
            "-d",
            "test",
            "-s",
            str(bad_scorefile),
            "-t",
            str(good_variants),
            "--outdir",
            str(outdir),
            "--min_overlap",
            "0.75",
        )
    ]
    flargs = list(itertools.chain(*args))

    with pytest.raises(ZeroMatchesError):
        with patch("sys.argv", flargs):
            run_match()
