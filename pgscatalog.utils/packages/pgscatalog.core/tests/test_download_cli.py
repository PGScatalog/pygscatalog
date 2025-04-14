import os
import pytest

from pgscatalog.core.cli.download_cli import run
from unittest.mock import patch

trait = "EFO_0004214"
pgs_id = "PGS000822"
pgp_id = "PGP000001"


def test_score_download(tmp_path):
    """Test downloading a scoring file with a PGS ID"""
    out_dir = str(tmp_path.resolve())

    args: list[str] = ["pgscatalog-download", "-i", pgs_id, "-o", out_dir, "-v"]

    with patch("sys.argv", args):
        run()
        hm_score_filename = f"{pgs_id}.txt.gz"
        assert hm_score_filename in os.listdir(out_dir)


@pytest.mark.parametrize("build", [("GRCh37"), ("GRCh38")])
def test_score_download_harmonised(tmp_path, build):
    """Test downloading harmonised scoring files in different genome builds"""
    out_dir = str(tmp_path.resolve())

    args: list[str] = [
        "pgscatalog-download",
        "-i",
        pgs_id,
        "--build",
        build,
        "-o",
        out_dir,
        "-v",
    ]

    with patch("sys.argv", args):
        run()
        hm_score_filename = f"{pgs_id}_hmPOS_{build}.txt.gz"
        assert hm_score_filename in os.listdir(out_dir)


def test_download_pub(tmp_path):
    """Test downloading scoring files with a PGP ID"""
    out_dir = str(tmp_path.resolve())

    args: list[str] = [
        "pgscatalog-download",
        "--pgp",
        pgp_id,
        "--build",
        "GRCh38",
        "-o",
        out_dir,
        "-v",
    ]

    with patch("sys.argv", args):
        run()
        assert len(os.listdir(out_dir)) == 3


def test_download_trait(tmp_path):
    """Test downloading scoring files with a trait"""
    out_dir = str(tmp_path.resolve())

    args: list[str] = [
        "pgscatalog-download",
        "--efo",
        trait,
        "--build",
        "GRCh38",
        "-o",
        out_dir,
        "-v",
    ]

    with patch("sys.argv", args):
        run()
        # 6 PGS for Abdominal Aortic Aneurysm 2024-01-16
        assert len(os.listdir(out_dir)) >= 6


def test_download_mix(tmp_path):
    """Test downloading scoring files with a mix of traits, publication, and score accessions"""
    out_dir = str(tmp_path.resolve())

    args: list[str] = [
        "pgscatalog-download",
        "--efo",
        trait,
        "--pgp",
        pgp_id,
        "-i",
        pgs_id,
        "--build",
        "GRCh38",
        "-o",
        out_dir,
        "-v",
    ]

    with patch("sys.argv", args):
        run()
        n_pgp = 3
        n_pgs = 1
        n_efo = 6
        # 6 PGS for Abdominal Aortic Aneurysm 2024-01-16
        assert len(os.listdir(out_dir)) >= n_pgp + n_pgs + n_efo
