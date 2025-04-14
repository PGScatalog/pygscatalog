import csv
import itertools
import os
from unittest.mock import patch

import pytest
from xopen import xopen

from pgscatalog.core.cli.relabel_cli import run


@pytest.fixture(scope="package")
def map_files(request):
    return [
        str(request.path.parent / "data" / "relabel_map_chr1.txt.gz"),
        str(request.path.parent / "data" / "relabel_map_chr8.txt.gz"),
    ]


@pytest.fixture(scope="package")
def relabel_scorefile(request):
    return str(request.path.parent / "data" / "hgdp_ALL_additive_0.scorefile")


def test_relabel(tmp_path_factory, relabel_scorefile, map_files):
    """Test relabelling a scorefile."""
    out_dir = tmp_path_factory.mktemp("outdir")

    args = [
        ("pgscatalog-relabel", "-m"),
        map_files,
        (
            "--target_file",
            relabel_scorefile,
            "--target_col",
            "ID",
            "-d",
            "hgdp",
            "--col_from",
            "ID_TARGET",
            "--col_to",
            "ID_REF",
            "-o",
            str(out_dir),
            "--combined",
            "--split",
        ),
    ]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run()

    out_f = os.listdir(out_dir)
    assert sorted(out_f) == [
        "hgdp_1_relabelled.gz",
        "hgdp_8_relabelled.gz",
        "hgdp_ALL_relabelled.gz",
    ]

    for x in out_f:
        with xopen(out_dir / x) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for line in reader:
                assert "ID" in line
                assert "effect_allele" in line
                assert "PGS000802_hmPOS_GRCh38" in line
