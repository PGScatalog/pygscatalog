import glob
import itertools
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
    return ["sampleset", "IID", "DENOM", "PGS", "SUM", "AVG"]


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
