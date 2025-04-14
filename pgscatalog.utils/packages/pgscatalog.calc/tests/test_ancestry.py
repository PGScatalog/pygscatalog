import gzip
import json
import os
from unittest.mock import patch

import pytest
import pandas as pd

from pgscatalog.calc.cli.ancestry_cli import run_ancestry


@pytest.fixture(scope="module")
def scorefile(request):
    return request.path.parent / "data" / "aggregated_scores.txt"


@pytest.fixture(scope="module")
def ref_psam(request):
    return request.path.parent / "data" / "ref.psam"


@pytest.fixture(scope="module")
def ref_pcs(request):
    return request.path.parent / "data" / "ref.pcs"


@pytest.fixture(scope="module")
def target_pcs(request):
    return request.path.parent / "data" / "target.pcs"


@pytest.fixture(scope="module")
def relatedness(request):
    return request.path.parent / "data" / "ref.king.cutoff.id"


def test_ancestry(
    scorefile, ref_psam, ref_pcs, target_pcs, relatedness, tmp_path_factory
):
    outdir = tmp_path_factory.mktemp("outdir")

    args = [
        "pgscatalog-ancestry-adjust",
        "-r",
        "reference",
        "--psam",
        str(ref_psam),
        "-s",
        str(scorefile),
        "--ref_pcs",
        str(ref_pcs),
        "--target_pcs",
        str(target_pcs),
        "-x",
        str(relatedness),
        "-p",
        "SuperPop",
        "--n_popcomp",
        "5",
        "-n",
        "empirical",
        "mean",
        "mean+var",
        "--n_normalization",
        "4",
        "--outdir",
        str(outdir),
        "--dataset",
        "hgdp",
        "--reference",
        "reference",
    ]

    with patch("sys.argv", args):
        run_ancestry()

    assert set(os.listdir(outdir)) == {
        "hgdp_info.json.gz",
        "hgdp_pgs.txt.gz",
        "hgdp_popsimilarity.txt.gz",
    }

    df = pd.read_table(outdir / "hgdp_pgs.txt.gz")
    assert set(df.columns) == {
        "Z_norm2",
        "percentile_MostSimilarPop",
        "Z_MostSimilarPop",
        "sampleset",
        "Z_norm1",
        "FID",
        "IID",
        "PGS",
        "SUM",
    }

    popsim = pd.read_table(outdir / "hgdp_popsimilarity.txt.gz")
    assert set(popsim.columns) == {
        "RF_P_EAS",
        "PC6",
        "FID",
        "IID",
        "PC4",
        "Unrelated",
        "PC2",
        "PC3",
        "PC1",
        "RF_P_AFR",
        "RF_P_SAS",
        "RF_P_EUR",
        "RF_P_AMR",
        "MostSimilarPop_LowConfidence",
        "PC10",
        "PC7",
        "MostSimilarPop",
        "sampleset",
        "SuperPop",
        "REFERENCE",
        "PC8",
        "PC9",
        "PC5",
    }

    with gzip.open(outdir / "hgdp_info.json.gz") as f:
        info = json.load(f)
        assert set(info.keys()) == {"pgs", "compare_pcs"}
