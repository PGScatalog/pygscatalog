import glob
import itertools
import shutil
from os import unlink
from unittest.mock import patch

from pgscatalog.calc.cli.aggregate_cli import run_aggregate

import pytest
import pandas as pd


@pytest.fixture(scope="module")
def scorefiles(request):
    basedir = request.path.parent / "data"
    return glob.glob(str(basedir / "hgdp*.zst"))


@pytest.fixture(scope="module")
def expected_output_columns():
    return ["sampleset", "FID", "IID", "PGS", "SUM", "DENOM", "AVG"]


def test_split_aggregate(tmp_path_factory, scorefiles, expected_output_columns):
    """Test aggregating HGDP PGS01229"""
    outdir = tmp_path_factory.mktemp("outdir")

    args = [
        ("pgscatalog-aggregate", "-s", *scorefiles, "--outdir", str(outdir), "--split")
    ]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run_aggregate()

    outf = list(outdir.glob("*.txt.gz"))
    assert [x.name for x in outf] == ["hgdp_pgs.txt.gz"]
    outdf = pd.read_csv(outf[0], sep="\t")
    assert list(outdf.columns) == expected_output_columns
    assert outdf.shape == (929, len(expected_output_columns))


def test_nosplit_aggregate(tmp_path_factory, scorefiles, expected_output_columns):
    """Test aggregating HGDP PGS01229 without splitting per sampleset"""
    outdir = tmp_path_factory.mktemp("outdir")

    args = [
        (
            "pgscatalog-aggregate",
            "-s",
            *scorefiles,
            "--outdir",
            str(outdir),
            "--no-split",
        )
    ]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run_aggregate()

    outf = list(outdir.glob("*.txt.gz"))
    assert [x.name for x in outf] == ["aggregated_scores.txt.gz"]
    outdf = pd.read_csv(outf[0], sep="\t")
    assert list(outdf.columns) == expected_output_columns
    assert outdf.shape == (929, len(expected_output_columns))


@pytest.fixture(scope="module")
def calculated_scores(request):
    """Calculated scores (split by effect type)"""
    basedir = request.path.parent / "data" / "aggregate"
    return glob.glob(str(basedir / "*sscore.zst"))


@pytest.fixture(scope="module")
def score_files(request):
    """Score variable reports from plink2"""
    basedir = request.path.parent / "data" / "aggregate"
    return glob.glob(str(basedir / "*scorefile.gz"))


@pytest.fixture(scope="module")
def score_vars(request):
    """Scoring files produced by pgscatalog-combine"""
    basedir = request.path.parent / "data" / "aggregate"
    return glob.glob(str(basedir / "*sscore.vars"))


def test_var_overlap(tmp_path_factory, calculated_scores, score_vars, score_files):
    outdir = tmp_path_factory.mktemp("outdir")

    args = [
        (
            "pgscatalog-aggregate",
            "-s",
            *calculated_scores,
            "--outdir",
            str(outdir),
            "--no-split",
            "--verify_variants",
        )
    ]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run_aggregate()

    outf = list(outdir.glob("*.txt.gz"))
    assert [x.name for x in outf] == ["aggregated_scores.txt.gz"]


@pytest.mark.parametrize(
    "use_score_vars, use_score_files",
    [
        (True, False),  # Use score_vars but not score_files
        (False, True),  # Use score_files but not score_vars
    ],
)
def test_var_overlap_missing(
    tmp_path_factory,
    calculated_scores,
    score_vars,
    score_files,
    use_score_files,
    use_score_vars,
):
    outdir = tmp_path_factory.mktemp("outdir")
    # stage data to outdir
    copied_calculated_scores = [shutil.copy(x, outdir) for x in calculated_scores]
    copied_score_vars = [shutil.copy(x, outdir) for x in score_vars]
    copied_score_files = [shutil.copy(x, outdir) for x in score_files]

    args = [
        (
            "pgscatalog-aggregate",
            "-s",
            *copied_calculated_scores,
            "--outdir",
            str(outdir),
            "--no-split",
            "--verify_variants",
            "--verbose",
        )
    ]
    flargs = list(itertools.chain(*args))

    # parameterised: delete data important for verification
    if not use_score_vars:
        [unlink(x) for x in copied_score_vars]

    if not use_score_files:
        [unlink(x) for x in copied_score_files]

    with pytest.raises(FileNotFoundError):
        with patch("sys.argv", flargs):
            run_aggregate()


def test_var_overlap_fails(
    tmp_path_factory, calculated_scores, score_vars, score_files
):
    outdir = tmp_path_factory.mktemp("outdir")
    # stage data to outdir
    copied_calculated_scores = [shutil.copy(x, outdir) for x in calculated_scores]
    copied_score_vars = [shutil.copy(x, outdir) for x in score_vars]
    _ = [shutil.copy(x, outdir) for x in score_files]

    with open(copied_score_vars[0], mode="at") as f:
        f.writelines(["potato"])

    args = [
        (
            "pgscatalog-aggregate",
            "-s",
            *copied_calculated_scores,
            "--outdir",
            str(outdir),
            "--no-split",
            "--verify_variants",
            "--verbose",
        )
    ]
    flargs = list(itertools.chain(*args))

    with pytest.raises(ValueError) as excinfo:
        with patch("sys.argv", flargs):
            run_aggregate()

    assert "Missing variants" in str(excinfo.value)
