import itertools
import os
from unittest.mock import patch

from pgscatalog.match.cli.intersect_cli import run_intersect
import pytest


@pytest.fixture(scope="module")
def vmiss(request):
    return request.path.parent / "data" / "target.vmiss.gz"


@pytest.fixture(scope="module")
def afreq(request):
    return request.path.parent / "data" / "target.afreq.gz"


@pytest.fixture(scope="module")
def target(request):
    return request.path.parent / "data" / "target.pvar.zst"


@pytest.fixture(scope="module")
def ref(request):
    return request.path.parent / "data" / "ref.pvar.zst"


def test_intersect_cli(tmp_path_factory, ref, target, afreq, vmiss):
    outdir = tmp_path_factory.mktemp("outdir")

    args = [
        (
            "pgscatalog-intersect",
            "-r",
            str(ref),
            "-t",
            str(target),
            "--outdir",
            str(outdir),
        )
    ]
    flargs = list(itertools.chain(*args))

    with patch("sys.argv", flargs):
        run_intersect()

    assert sorted(os.listdir(outdir)) == [
        "intersect_counts_None.txt",
        "matched_variants.txt",
        "reference_variants.txt",
        "target_variants.txt",
    ]

    # header + 100 variants
    with open(outdir / "target_variants.txt") as f:
        assert sum(1 for _ in f) == 101

    with open(outdir / "reference_variants.txt") as f:
        # included 100 extra variants in reference that won't intersect intentionally
        assert sum(1 for _ in f) == 201

    with open(outdir / "matched_variants.txt") as f:
        assert sum(1 for _ in f) == 101

    with open(outdir / "intersect_counts_None.txt") as f:
        intersect_counts = f.read().splitlines()
        # 100 variants in target
        # 200 variants in reference
        # 100 variants intersected successfully
        assert intersect_counts == ["100", "200", "100"]
