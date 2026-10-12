"""A serial run's peak memory is the largest single stage, not the sum of four.

This is the claim that makes `pgsc_calc run` usable on real data: SCORE still
holds one shard (ADR 0003) and MATCH still sees one chromosome, so what a serial
run gives up is wall-clock and nothing else. It holds only if each stage releases
what it held - Nextflow gets that for free by ending the process, and one leaked
store, open duckdb connection or retained dask graph turns max into sum.

Every measurement runs in its own interpreter. `ru_maxrss` is a high-water mark
that never falls, so read inside the pytest process it would report whatever the
rest of the suite peaked at rather than what a stage cost.
"""

import pathlib
import shutil
import subprocess
import sys

from pgscatalog.calc.cli import run

from .conftest import SAMPLESET, run_cli, target_spec

_MEASURE = """
import sys, resource
from pgscatalog.calc.cli.main import main

main(sys.argv[1:])
print("peak_rss", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
"""


def _peak_rss(argv: list[str]) -> int:
    """Run one command in its own interpreter and report what it peaked at

    Units are bytes on macOS and kilobytes on Linux. Both sides of the comparison
    come from the same platform, so the ratio cancels them - don't "fix" this by
    scaling one of them.
    """
    result = subprocess.run(
        [sys.executable, "-c", _MEASURE, *argv],
        check=True,
        capture_output=True,
        text=True,
    )
    peaks = [x for x in result.stdout.splitlines() if x.startswith("peak_rss ")]
    return int(peaks[-1].split()[1])


def _stage_dir(argv: list[str]) -> str:
    """Where a recorded stage wrote, which has to be emptied before it reruns"""
    flag = "--cache_dir" if argv[0] == "load" else "--out_dir"
    return argv[argv.index(flag) + 1]


def test_a_serial_run_peaks_at_its_largest_stage(
    grid_targets, scorefile, tmp_path, monkeypatch
):
    """Every stage releasing what it held is what makes a serial run's limit the
    number of shards rather than the size of the dataset (ADR 0013).
    """

    def serial_argv(out_dir) -> list[str]:
        return [
            "run",
            "--target_genomes",
            *[target_spec(x) for x in grid_targets],
            "--sampleset",
            SAMPLESET,
            "--format",
            "vcf",
            "--score_paths",
            str(scorefile),
            "--out_dir",
            str(out_dir),
        ]

    # `run` composes the argv of every stage it executes, so recording what it
    # dispatched measures each stage on exactly the arguments it ran with, rather
    # than rebuilding the stage graph here a second time
    stages: list[list[str]] = []
    dispatch = run._dispatch

    def record(stage_argv: list[str], dispatched: list[str]) -> None:
        stages.append(list(stage_argv))
        dispatch(stage_argv, dispatched)

    monkeypatch.setattr(run, "_dispatch", record)
    run_cli(["pgsc_calc", *serial_argv(tmp_path / "recorded")], monkeypatch)
    monkeypatch.undo()

    # each stage reruns against the tree that run left behind, into its own
    # directory emptied first, so it does exactly the same work in a clean process
    peaks = []
    for stage_argv in stages:
        stage_dir = _stage_dir(stage_argv)
        shutil.rmtree(stage_dir)
        peaks.append(
            (f"{stage_argv[0]} {pathlib.Path(stage_dir).name}", _peak_rss(stage_argv))
        )

    serial = _peak_rss(serial_argv(tmp_path / "measured"))
    largest = max(peak for _, peak in peaks)

    # measured on grid_targets (macOS, 4 shards): LOAD 397 MB, MATCH 303 MB,
    # MATCH_COMBINE 299 MB, SCORE 432 MB, SCORE_COMBINE 274 MB, and the whole
    # serial run 551 MB - 1.3x its largest stage, against 1.7 GB if it summed.
    # Every measurement also carries a fresh interpreter's baseline, which is on
    # both sides and so dilutes the ratio rather than biasing it. The factor is
    # slack around a regression into `sum`, not a number worth pinning.
    assert serial <= 2 * largest, (
        f"a serial run peaked at {serial}, its largest stage at {largest}: {peaks}"
    )
