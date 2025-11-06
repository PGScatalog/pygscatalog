from __future__ import annotations

import argparse
import atexit
import logging
import logging.handlers
import pathlib
from dataclasses import dataclass
from shutil import rmtree

from rich.logging import RichHandler
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TaskID,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

from pgscatalog.calc import (
    VALID_CHROMOSOMES,
    GenomeFileType,
    Pathish,
    Scorefiles,
    TargetGenome,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class AnnotatedGenomePath:
    path: str  # intentionally a str to support s3
    chrom: str | None

    @classmethod
    def from_string(cls, s: str) -> AnnotatedGenomePath:
        try:
            path, chrom = s.rsplit(":", 1)
        except ValueError:
            path, chrom = s, None

        return cls(path=path, chrom=chrom)


def parse_genome_paths(
    target_genomes: list[str],
    cache_dir: Pathish,
    sampleset: str,
    bgen_sample_file: Pathish | None,
) -> list[TargetGenome]:
    targets: list[AnnotatedGenomePath] = list(
        {AnnotatedGenomePath.from_string(x) for x in target_genomes}
    )

    genomes = []

    for genome in targets:
        genomes.append(
            TargetGenome(
                target_path=genome.path,
                chrom=genome.chrom,
                cache_dir=cache_dir,
                sampleset=sampleset,
                sample_file=bgen_sample_file,
            )
        )

    return genomes


def check_chromosomes(target_genomes: list[TargetGenome]) -> None:
    if not {x.chrom for x in target_genomes}.issubset(VALID_CHROMOSOMES):
        raise NotImplementedError("Chromosomes 1-22 supported but this will be fixed")

    if all(x.chrom is not None for x in target_genomes):
        is_chrom_split = True
    else:
        is_chrom_split = False

    logger.info(f"{is_chrom_split=}")

    if not is_chrom_split and len(target_genomes) > 1:
        logger.critical(
            "Unsplit chromosome detected but more than one VCF is set in "
            f"--target_genomes {target_genomes=}"
        )
        raise TypeError(
            "Missing chromosome information in file path (try --help for details)"
        )


def load_cli(args: argparse.Namespace) -> None:
    if args.format == GenomeFileType.BGEN and not args.bgen_sample_file:
        raise argparse.ArgumentTypeError(
            "--bgen_sample_file is required when --format is 'bgen'"
        )

    args.cache_dir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        level="WARNING",
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True)],
    )

    if args.verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

    atexit.register(cleanup_cache_tmpdir, cache_dir=args.cache_dir)

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        MofNCompleteColumn(),
        BarColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        refresh_per_second=2,
    ) as progress:
        progress.print(":test_tube:", "pgsc_calc load", ":dna:", "starting...")

        # set up input data
        target_genomes: list[TargetGenome] = parse_genome_paths(
            target_genomes=args.target_genomes,
            cache_dir=args.cache_dir,
            sampleset=args.sampleset,
            bgen_sample_file=args.bgen_sample_file,
        )

        # check
        _is_split: bool = check_chromosomes(target_genomes=target_genomes)
        scorefiles: Scorefiles = Scorefiles(args.score_paths)

        for target_genome in target_genomes:
            positions_to_query = scorefiles.get_unique_positions(
                chrom=target_genome.chrom, zarr_group=target_genome.zarr_group
            )

            if len(positions_to_query) == 0:
                progress.print("All variants are already cached. Yay!")
                continue
            target_genome.cache_variants(positions_to_query)

        progress.print("Finished caching :tada: Goodbye!")


def run_jobs_sequentially(
    progress: Progress,
    task_id: TaskID,
    jobs: list[tuple[TargetGenome, tuple[tuple[str, int], ...]]],
) -> None:
    progress.print(f"Running {len(jobs)} jobs sequentially")
    for i, job in enumerate(jobs, start=1):
        target_genome, positions = job
        progress.print(
            f"Job {i}/{len(jobs)} started: {len(positions)} positions to cache"
        )
        target_genome.cache_variants(positions=positions)
        progress.update(task_id, completed=i, refresh=True)


def cleanup_cache_tmpdir(cache_dir: pathlib.Path) -> None:
    # use print statements in this atexit handler
    # during interpreter exit logging handlers might not exist
    tmp_dir = cache_dir / "tmp"
    print(f"Cleaning up {tmp_dir} if it exists")
    try:
        rmtree(tmp_dir)
    except FileNotFoundError:
        print("Nothing to clean up")
    finally:
        print("cleanup function exiting, goodbye!")
