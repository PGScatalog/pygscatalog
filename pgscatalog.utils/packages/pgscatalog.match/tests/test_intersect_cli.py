import functools
import itertools
import os
import gzip
import shutil
import tempfile
from unittest.mock import patch

from pgscatalog.match.cli.intersect_cli import read_var_general, run_intersect
import pytest


@pytest.fixture(scope="module")
def vmiss(request):
    return request.path.parent / "data" / "target.vmiss.gz"


@pytest.fixture(scope="module")
def afreq(request):
    return request.path.parent / "data" / "target.afreq.gz"


@pytest.fixture(scope="module")
def ref(request):
    return request.path.parent / "data" / "ref.pvar.zst"


@pytest.fixture(params=["target.bim.zst", "target.pvar.zst"], scope="module")
def target(request):
    return request.path.parent / "data" / request.param


def test_intersect_cli(tmp_path_factory, ref, target, afreq, vmiss):
    outdir = tmp_path_factory.mktemp("outdir")

    # using a low batch size to test batching behaviour works as expected
    args = [
        (
            "pgscatalog-intersect",
            "-r",
            str(ref),
            "-t",
            str(target),
            "--outdir",
            str(outdir),
            "--batch_size",
            "100",
        )
    ]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run_intersect()

    assert sorted(os.listdir(outdir)) == [
        "intersect_counts_None.txt",
        "matched_variants.txt.gz",
        "reference_variants.txt.gz",
        "target_variants.txt.gz",
    ]

    # header + 100 variants
    with gzip.open(outdir / "target_variants.txt.gz") as f:
        assert sum(1 for _ in f) == 101

    with gzip.open(outdir / "reference_variants.txt.gz") as f:
        # included 100 extra variants in reference that won't intersect intentionally
        assert sum(1 for _ in f) == 201

    with gzip.open(outdir / "matched_variants.txt.gz") as f:
        assert sum(1 for _ in f) == 101

    with open(outdir / "intersect_counts_None.txt") as f:
        intersect_counts = f.read().splitlines()
        # 100 variants in target
        # 200 variants in reference
        # 100 variants intersected successfully
        assert intersect_counts == ["100", "200", "100"]


def count_lines(path):
    with gzip.open(path) as f:
        return sum(1 for _ in f)


def test_intersect_cli_sort_chunk_name_contains_bim(tmp_path, ref, target):
    # sort chunks are random temporary files, which are headered tsv and not bim
    args = ["pgscatalog-intersect", "-r", str(ref), "-t", str(target)]
    args += ["--outdir", str(tmp_path), "--batch_size", "100"]
    chunk = functools.partial(tempfile.NamedTemporaryFile, prefix="tmpbim")

    with patch("sys.argv", args), patch("tempfile.NamedTemporaryFile", chunk):
        run_intersect()

    assert count_lines(tmp_path / "matched_variants.txt.gz") == 101


def test_intersect_cli_pvar_name_contains_bim(tmp_path, ref):
    for extension in ("pvar.zst", "afreq.gz", "vmiss.gz"):
        shutil.copy(
            ref.parent / f"target.{extension}", tmp_path / f"testbim_22.{extension}"
        )
    target = tmp_path / "testbim_22.pvar.zst"
    args = ["pgscatalog-intersect", "-r", str(ref), "-t", str(target)]
    args += ["--outdir", str(tmp_path)]

    with patch("sys.argv", args):
        run_intersect()

    assert count_lines(tmp_path / "matched_variants.txt.gz") == 101


def test_read_var_general_uncompressed_bim(tmp_path):
    bim = tmp_path / "target.bim"
    bim.write_text("1\trs1\t0\t100\tA\tG\n")

    assert next(read_var_general(bim))["ID"] == "rs1"
