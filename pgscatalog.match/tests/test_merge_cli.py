import itertools
from unittest.mock import patch
import pytest

from glob import glob

from pgscatalog.core.lib.pgsexceptions import ZeroMatchesError
from pgscatalog.match.cli.merge_cli import run_merge


@pytest.fixture(scope="module")
def good_scorefile(request):
    return request.path.parent / "data" / "good_match_scorefile.txt"


@pytest.fixture(scope="module")
def bad_scorefile(request):
    return request.path.parent / "data" / "bad_match_scorefile.txt.gz"


@pytest.fixture(scope="module")
def match_ipc(request):
    return request.path.parent / "data" / "match.ipc"


def test_merge(tmp_path_factory, good_scorefile, match_ipc):
    """Test merging runs without errors with good data"""
    outdir = tmp_path_factory.mktemp("outdir")

    args = [
        (
            "pgscatalog-matchmerge",
            "-d",
            "test",
            "-s",
            str(good_scorefile),
            "--matches",
            str(match_ipc),
            "--outdir",
            str(outdir),
            "--min_overlap",
            str(0.75),
        )
    ]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run_merge()

    assert (outdir / "test_log.csv.gz").exists()
    assert (outdir / "test_summary.csv").exists()
    assert (outdir / "test_ALL_recessive_0.scorefile.gz").exists()
    assert (outdir / "test_ALL_dominant_0.scorefile.gz").exists()
    assert (outdir / "test_ALL_additive_0.scorefile.gz").exists()


def test_split_output(tmp_path_factory, good_scorefile, match_ipc):
    """Test merging runs without errors with good data outputting split data"""
    outdir = tmp_path_factory.mktemp("outdir")

    args = [
        (
            "pgscatalog-matchmerge",
            "-d",
            "test",
            "-s",
            str(good_scorefile),
            "--matches",
            str(match_ipc),
            "--outdir",
            str(outdir),
            "--min_overlap",
            str(0.75),
            "--split",
        )
    ]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run_merge()

    # 1 dominant, 1 recessive, 19 additive
    splitf = glob(str(outdir / "*scorefile.gz"))
    assert len(splitf) == 21
    assert sum("additive" in x for x in splitf) == 19
    assert sum("dominant" in x for x in splitf) == 1
    assert sum("recessive" in x for x in splitf) == 1


def test_splitcombine_output(tmp_path_factory, good_scorefile, match_ipc):
    """Test merging runs without errors with good data outputting both combined and split scorefiles"""
    outdir = tmp_path_factory.mktemp("outdir")

    args = [
        (
            "pgscatalog-matchmerge",
            "-d",
            "test",
            "-s",
            str(good_scorefile),
            "--matches",
            str(match_ipc),
            "--outdir",
            str(outdir),
            "--min_overlap",
            str(0.75),
            "--split",
            "--combined",
        )
    ]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run_merge()
    splitf = glob(str(outdir / "*scorefile.gz"))
    assert len(splitf) == 21 + 3  # split + combined
    assert sum("additive" in x for x in splitf) == 19 + 1
    assert sum("dominant" in x for x in splitf) == 1 + 1
    assert sum("recessive" in x for x in splitf) == 1 + 1


def test_strict_merge(tmp_path_factory, good_scorefile, match_ipc):
    """Test merging with extremely strict overlap to trigger a ZeroMatchesError"""
    outdir = tmp_path_factory.mktemp("outdir")

    args = [
        (
            "pgscatalog-matchmerge",
            "-d",
            "test",
            "-s",
            str(good_scorefile),
            "--matches",
            str(match_ipc),
            "--outdir",
            str(outdir),
            "--min_overlap",
            str(1.00),
        )
    ]
    flargs = list(itertools.chain(*args))

    with pytest.raises(ZeroMatchesError):
        with patch("sys.argv", flargs):
            run_merge()

    # don't write any scoring files
    assert glob(str(outdir / "*scorefile.gz")) == []
