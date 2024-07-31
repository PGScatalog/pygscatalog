import csv
import gzip
import itertools
import os
import pathlib
from unittest.mock import patch
import pytest

from pgscatalog.core import ZeroMatchesError
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


@pytest.fixture(scope="module")
def multiallelic_variants(request):
    return request.path.parent / "data" / "multiallelic.pvar"


@pytest.fixture()
def filter_ids(request, tmp_path, good_scorefile):
    with open(good_scorefile) as f:
        scorefile = list(csv.DictReader(f, delimiter="\t"))
    # build variant IDs from the scoring file
    ids = [
        ":".join(
            [x["chr_name"], x["chr_position"], x["effect_allele"], x["other_allele"]]
        )
        for x in scorefile
    ]
    # grab half of them, use these to filter matches
    filter_list = ids[: int(len(ids) / 2)]
    outf = tmp_path / "filter_ids.txt"
    with open(outf, mode="wt") as f:
        [f.write(x + "\n") for x in filter_list]

    return pathlib.Path(outf)


def test_filter_match(tmp_path_factory, good_variants, good_scorefile, filter_ids):
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
            "--filter_IDs",
            str(filter_ids),
        )
    ]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run_match()

    with gzip.open(
        outdir / "test_ALL_additive_0.scorefile.gz", mode="rt"
    ) as out_f, open(good_scorefile) as f1, open(filter_ids) as f2:
        n_out_score = len(list(csv.DictReader(out_f, delimiter="\t")))
        n_scorefile = sum(1 for _ in f1)
        n_filt = sum(1 for _ in f2)

    # matched output from scorefile has been filtered
    assert n_out_score < n_filt < n_scorefile


def test_multiallelic(tmp_path_factory, multiallelic_variants, good_scorefile):
    outdir = tmp_path_factory.mktemp("outdir")

    args = [
        (
            "pgscatalog-match",
            "-d",
            "test",
            "-s",
            str(good_scorefile),
            "-t",
            str(multiallelic_variants),
            "--outdir",
            str(outdir),
            "--keep_multiallelic",
            "--min_overlap",
            "0.75",
        )
    ]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run_match()

    assert (outdir / "test_summary.csv").exists()

    # test multiallelic matches happened
    with open(outdir / "test_summary.csv") as f:
        reader = csv.DictReader(f)
        log = list(reader)

    assert any({"is_multiallelic": "true"}.items() <= x.items() for x in log)
    assert (outdir / "test_log.csv.gz").exists()
    assert (outdir / "test_ALL_recessive_0.scorefile.gz").exists()
    assert (outdir / "test_ALL_dominant_0.scorefile.gz").exists()
    assert (outdir / "test_ALL_additive_0.scorefile.gz").exists()


@pytest.mark.parametrize(
    "ambiguous,multiallelic,skipflip,keep_first_match",
    [
        (
            "--keep_ambiguous",
            "--keep_multiallelic",
            "--ignore_strand_flips",
            "--keep_first_match",
        ),
        (None, None, None, None),
    ],
)
def test_match(
    tmp_path_factory,
    good_scorefile,
    good_variants,
    ambiguous,
    multiallelic,
    skipflip,
    keep_first_match,
):
    """Test matching runs without errors with good data

    Parameterised to run twice: with default CLI match args and optional matching arguments"""
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
            ambiguous,
            multiallelic,
            skipflip,
            keep_first_match,
        )
    ]
    flargs = list(itertools.chain(*args))
    flargs = [x for x in flargs if x]  # drop parameterised None

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


def test_stringent_match(tmp_path_factory, good_scorefile, good_variants):
    """Test matching well with stringent match rate.

    Some scores pass, some fail, but the application exits successfully.
    """
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
            "0.9",
        )
    ]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run_match()

    assert (outdir / "test_log.csv.gz").exists()
    assert (outdir / "test_summary.csv").exists()

    with open(outdir / "test_summary.csv") as f:
        summary = list(csv.DictReader(f))
        # at least one score fails matching
        assert any([x["score_pass"] == "false" for x in summary])
        # and at least one score passes matching
        assert any([x["score_pass"] == "true" for x in summary])

    # but scoring files are still written because at least one score passed
    assert (outdir / "test_ALL_additive_0.scorefile.gz").exists()


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
