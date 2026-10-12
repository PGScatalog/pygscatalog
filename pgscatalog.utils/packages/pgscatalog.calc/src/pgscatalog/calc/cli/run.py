"""A Serial run: LOAD every target, then the four scoring stages in sequence.

Nextflow dispatches those stages as separate tasks (ADR 0001). This runs the same
stages, in the same order, over the same on-disk artefacts, with nothing
dispatching them - the entry point for environments where Nextflow can't run
(ADR 0013). It is a convenience, not a second scaling mechanism: SCORE still holds
one shard and MATCH still sees one chromosome, so peak memory is the largest
single stage and what a serial run gives up is wall-clock. Its limit is therefore
the number of shards, not the size of the dataset.
"""

from __future__ import annotations

import logging
import shlex
import shutil
from typing import TYPE_CHECKING

from rich import print

from pgscatalog.calc.cli._stage import add_stage_args, configure_logging
from pgscatalog.calc.cli._utils import zero_one_float
from pgscatalog.calc.cli.load import (
    AnnotatedGenomePath,
    add_load_input_args,
    check_leaf_paths,
)
from pgscatalog.calc.lib.scorefile import Scorefiles

if TYPE_CHECKING:
    import argparse
    import pathlib

logger = logging.getLogger(__name__)

STAGES = ("load", "match", "match-combine", "score", "combine")
"""Every stage a serial run executes.

Declared, and checked against what `run` actually dispatched, so that dropping a
stage is a choice rather than an oversight. Nothing here can detect a *fifth*
stage added to main.nf and never wired in - that risk is accepted (ADR 0013).
"""


def add_run_args(parser: argparse.ArgumentParser) -> None:
    add_load_input_args(parser)
    parser.add_argument(
        "--min_overlap",
        type=zero_one_float,
        dest="min_overlap",
        help="Minimum variant overlap",
        required=False,
        default=0.75,
    )
    add_stage_args(parser)
    parser.set_defaults(func=run_cli)


def _dispatch(argv: list[str], dispatched: list[str]) -> None:
    """Run one stage by parsing its argv, exactly as its own subcommand would

    Composing argv instead of hand-building a Namespace makes argparse the drift
    check: a renamed or newly-required stage argument fails naming the flag,
    rather than raising AttributeError part-way through a stage that has already
    done an hour of work (ADR 0013).
    """
    # imported here because main imports this module to build its subparsers
    from pgscatalog.calc.cli.main import build_parser

    command = shlex.join(["pgsc_calc", *argv])
    logger.info(f"Running stage: {command}")

    try:
        args = build_parser().parse_args(argv)
        args.func(args)
    except Exception:
        # there is no resume: the completed stages' outputs are left on disk and
        # this command reruns the one that failed against them
        logger.error(f"Stage failed. Rerun it with: {command}")
        raise

    dispatched.append(argv[0])


def run_cli(args: argparse.Namespace) -> None:
    configure_logging(args.verbose)

    out_dir: pathlib.Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    stages_dir = out_dir / "stages"

    dispatched: list[str] = []
    score_paths = [str(x) for x in args.score_paths]
    resources = [
        "--threads",
        str(args.threads),
        "--max_memory_gb",
        str(args.max_memory_gb),
    ]
    verbose = ["--verbose"] if args.verbose else []
    sample_file = (
        ["--bgen_sample_file", str(args.bgen_sample_file)]
        if args.bgen_sample_file
        else []
    )
    # without this `run` cannot score any X-carrying PGS, because LOAD refuses one
    # without a declaration (ADR 0014)
    sex_map = ["--sex_map", str(args.sex_map)] if args.sex_map else []

    print(":test_tube:", "pgsc_calc run", ":dna:", "starting...")

    # LOAD once per target, into its own cache: SCORE takes one zip holding one
    # shard, so one archive per target genome is a requirement (ADR 0003)
    zips: list[pathlib.Path] = []
    chromosomes: list[str] = []

    targets = [
        AnnotatedGenomePath.from_string(target) for target in args.target_genomes
    ]
    check_leaf_paths(targets, args.sampleset)
    score_chromosomes = {
        chrom for chrom, _ in Scorefiles(args.score_paths).get_unique_positions()
    }
    for i, (target, annotated) in enumerate(
        zip(args.target_genomes, targets, strict=True)
    ):
        if annotated.chrom is not None and annotated.chrom not in score_chromosomes:
            logger.info(f"Skipping {target}: no scoring positions on this chromosome")
            continue
        cache_dir = stages_dir / "load" / str(i)
        _dispatch(
            [
                "load",
                "--target_genomes",
                target,
                "--sampleset",
                args.sampleset,
                "--format",
                args.format.value,
                "--score_paths",
                *score_paths,
                "--absence_semantics",
                args.absence_semantics,
                "--batch_size",
                str(args.batch_size),
                "--cache_dir",
                str(cache_dir),
                "--threads",
                str(args.threads),
                *sample_file,
                *sex_map,
                *verbose,
            ],
            dispatched,
        )

        zip_path = cache_dir / "genotypes.zarr.zip"
        if not zip_path.exists():
            # LOAD writes no archive when nothing it queried was cached, which is
            # a shard the later stages would never hear about
            raise FileNotFoundError(
                f"{target} cached no variants, so {zip_path} doesn't exist. Check "
                "the scoring file positions and genome build match this target."
            )

        zips.append(zip_path)
        # the same parse LOAD uses, so chr1 and 1 group together exactly as they
        # land in one leaf of the cache
        chromosomes.append(annotated.chrom or "ALL")

    if not zips:
        raise ValueError("No target chromosomes have requested scoring positions")

    # MATCH sees one chromosome whole: every shard of it, and only it
    shards_by_chromosome: dict[str, list[pathlib.Path]] = {}
    for chromosome, zip_path in zip(chromosomes, zips):
        shards_by_chromosome.setdefault(chromosome, []).append(zip_path)

    matches: list[str] = []
    samples: list[str] = []

    for chromosome, chromosome_zips in shards_by_chromosome.items():
        match_dir = stages_dir / "match" / chromosome
        _dispatch(
            [
                "match",
                "--zarr_zip_file",
                *[str(x) for x in chromosome_zips],
                "--score_paths",
                *score_paths,
                "--batch_size",
                str(args.batch_size),
                "--out_dir",
                str(match_dir),
                *resources,
                *verbose,
            ],
            dispatched,
        )
        matches.append(str(match_dir / "matches.parquet"))
        samples.append(str(match_dir / "samples.parquet"))

    match_combine_dir = stages_dir / "match_combine"
    _dispatch(
        [
            "match-combine",
            "--matches",
            *matches,
            "--samples",
            *samples,
            "--score_paths",
            *score_paths,
            "--sampleset",
            args.sampleset,
            "--min_overlap",
            str(args.min_overlap),
            "--out_dir",
            str(match_combine_dir),
            *resources,
            *verbose,
        ],
        dispatched,
    )

    partial_scores: list[pathlib.Path] = []

    for i, zip_path in enumerate(zips):
        score_dir = stages_dir / "score" / str(i)
        _dispatch(
            [
                "score",
                "--zarr_zip_file",
                str(zip_path),
                "--score_paths",
                *score_paths,
                "--match_table",
                str(match_combine_dir / "match_table.parquet"),
                "--batch_size",
                str(args.batch_size),
                "--out_dir",
                str(score_dir),
                *resources,
                *verbose,
            ],
            dispatched,
        )
        partial_scores.append(score_dir / "partial_scores.parquet")

    # SCORE_COMBINE sums the partial scores it's handed and can't know how many
    # shards existed, so a missing one comes back as a quietly low score for every
    # sample. run is the only component that knows the expected count.
    if missing := [str(x) for x in partial_scores if not x.exists()]:
        raise FileNotFoundError(
            f"{len(partial_scores)} shards were scored but these partial scores "
            f"are missing, so the final scores would be too low: {missing}"
        )

    combine_dir = stages_dir / "combine"
    _dispatch(
        [
            "combine",
            "--partial_scores",
            *[str(x) for x in partial_scores],
            "--out_dir",
            str(combine_dir),
            *resources,
            *verbose,
        ],
        dispatched,
    )

    scores = combine_dir / "scores.txt.gz"
    shutil.copy2(scores, out_dir / scores.name)

    if set(dispatched) != set(STAGES):
        raise ValueError(
            f"A serial run must execute every stage in {STAGES}, but it ran "
            f"{sorted(set(dispatched))}"
        )

    print(f"Finished serial run :tada: Scores are in {out_dir / scores.name}")
