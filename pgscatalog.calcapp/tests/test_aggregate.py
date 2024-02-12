import itertools
from unittest.mock import patch

from pgscatalog.calcapp.aggregate_cli import run_aggregate

import pytest
import pandas as pd


@pytest.fixture(scope="module")
def scorefiles(request):
    return [str(x) for x in list(request.path.parent.glob("testdata/*.zst"))]


def test_split_aggregate(tmp_path_factory, scorefiles):
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
    assert list(outdf.columns) == [
        "sampleset",
        "IID",
        "DENOM",
        "PGS001229_hmPOS_GRCh38_SUM",
    ]
    assert outdf.shape == (929, 4)


def test_nosplit_aggregate(tmp_path_factory, scorefiles):
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
    assert list(outdf.columns) == [
        "sampleset",
        "IID",
        "DENOM",
        "PGS001229_hmPOS_GRCh38_SUM",
    ]
    assert outdf.shape == (929, 4)
