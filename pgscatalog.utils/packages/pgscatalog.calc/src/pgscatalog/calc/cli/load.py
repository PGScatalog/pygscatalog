from __future__ import annotations

import argparse
import atexit
import concurrent.futures
import functools
import itertools
import logging
import logging.handlers
import multiprocessing
import pathlib
from dataclasses import dataclass
from random import sample
from shutil import rmtree
from threading import Thread
from typing import TYPE_CHECKING

import duckdb
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

from pgscatalog.calc import GenomeFileType, Pathish, Scorefiles, TargetGenome

if TYPE_CHECKING:
    from collections.abc import Sequence
    from queue import Queue


logger = logging.getLogger("pgscatalog.wgs")


def logger_thread(queue: Queue[logging.LogRecord | None]) -> None:
    """This function runs in a thread in the main process"""
    while True:
        record: logging.LogRecord | None = queue.get()
        if record is None:
            break
        else:
            thread_log = logging.getLogger(record.name)
            thread_log.handle(record)


def init_worker(queue: Queue[logging.LogRecord | None], verbose: bool) -> None:
    """Send all worker logs to a queue.
    Worker processes must be initialised with this function.
    See https://docs.python.org/3/howto/logging-cookbook.html#logging-to-a-single-file-from-multiple-processes
    """
    h = logging.handlers.QueueHandler(queue)
    root = logging.getLogger()
    root.addHandler(h)
    if verbose:
        root.setLevel(logging.INFO)
    else:
        root.setLevel(logging.WARNING)


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
) -> dict[str | None, list[TargetGenome]]:
    targets: list[AnnotatedGenomePath] = list(
        set(AnnotatedGenomePath.from_string(x) for x in target_genomes)
    )

    genomes: dict[str | None, list[TargetGenome]] = {}

    for genome in targets:
        genomes.setdefault(genome.chrom, [])
        genomes[genome.chrom].append(
            TargetGenome(
                target_path=genome.path,
                chrom=genome.chrom,
                cache_dir=cache_dir,
                sampleset=sampleset,
                sample_file=bgen_sample_file,
            )
        )

    return genomes


def cache_variants(
    target: TargetGenome,
    positions: Sequence[tuple[str, int]],
) -> None:
    """A worker process function to fetch variants and store them in parquet files"""
    target.cache_variants(positions=positions)


def is_vcfs_split(target_genomes: dict[str | None, list[TargetGenome]]) -> bool:
    if all(k is not None for k, _ in target_genomes.items()):
        is_vcfs_split_by_chrom = True
    else:
        is_vcfs_split_by_chrom = False

    logger.info(f"{is_vcfs_split_by_chrom=}")

    if not is_vcfs_split_by_chrom and len(target_genomes) > 1:
        logger.critical(
            "Unsplit VCF chromosome detected but more than one VCF is set in "
            "--target_genomes"
        )
        raise TypeError(
            "Missing chromosome information in file path (try --help for details)"
        )

    return is_vcfs_split_by_chrom


def get_jobs(
    is_vcfs_split_by_chrom: bool,
    target_genomes: dict[str | None, list[TargetGenome]],
    unique_positions: list[tuple[str, int]],
    batch_size: int,
) -> list[tuple[TargetGenome, tuple[tuple[str, int], ...]]]:
    cache_jobs = []
    if is_vcfs_split_by_chrom:
        for k, g in itertools.groupby(
            sorted(unique_positions, key=lambda x: (x[0], x[1])), key=lambda p: p[0]
        ):
            genomes: list[TargetGenome] | None = target_genomes.get(k, None)
            if genomes is None:
                logger.warning(
                    f"Skipping scoring file positions in chromosome {k}: no associated"
                    f" target genome found"
                )
                continue
            else:
                variants = list(itertools.batched(g, batch_size))
                cache_jobs.extend(list(itertools.product(genomes, variants)))
    else:
        genomes = target_genomes.get(None, None)
        if genomes is None:
            raise TypeError("Can't get target genome")

        variants = list(itertools.batched(unique_positions, batch_size))
        cache_jobs.extend(list(itertools.product(genomes, variants)))

    # jobs are ordered by filename: shuffle them to avoid thundering herds
    return sample(cache_jobs, k=len(cache_jobs))


def get_cached_positions(cache_path: Pathish) -> set[tuple[str, int]]:
    # get cached positions from the DB before launching workers
    db_path = pathlib.Path(cache_path) / "variants.db"
    if db_path.exists():
        with duckdb.connect(str(db_path)) as db:
            pos = set(
                db.execute(
                    "SELECT DISTINCT chr_name, chr_pos FROM variants_table;"
                ).fetchall()
            )
    else:
        pos = set()
    return pos


def get_query_positions(
    scorefiles: Scorefiles, cache_dir: Pathish
) -> list[tuple[str, int]]:
    unique_positions: list[tuple[str, int]] = scorefiles.get_unique_positions()
    cached_positions = get_cached_positions(cache_path=cache_dir)
    return list(set(unique_positions) - cached_positions)


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

    # an unbounded queue suitable for concurrent.futures.ProcessPoolExecutor
    # see https://docs.python.org/3/howto/logging-cookbook.html#using-concurrent-futures-processpoolexecutor  # noqa: E501
    log_queue: Queue = multiprocessing.Manager().Queue(-1)
    lp = Thread(target=logger_thread, args=(log_queue,))
    lp.start()

    if args.verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        MofNCompleteColumn(),
        BarColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        refresh_per_second=2,
    ) as progress:
        progress.print(":test_tube:", "pgscatalog-genome-cache", ":dna:", "starting...")
        atexit.register(cleanup_cache_tmpdir, cache_dir=args.cache_dir)

        # set up input data
        target_genomes: dict[str | None, list[TargetGenome]] = parse_genome_paths(
            target_genomes=args.target_genomes,
            cache_dir=args.cache_dir,
            sampleset=args.sampleset,
            bgen_sample_file=args.bgen_sample_file,
        )
        is_vcfs_split_by_chrom: bool = is_vcfs_split(target_genomes=target_genomes)
        scorefiles: Scorefiles = Scorefiles(args.score_paths)
        positions_to_query = get_query_positions(
            scorefiles=scorefiles, cache_dir=args.cache_dir
        )

        if len(positions_to_query) == 0:
            progress.print("All variants are already cached. Yay!")
            return

        # set up caching jobs to parallelise work by creating batches of non-overlapping
        # positions. often creates multiple jobs for each file: helpful because
        # chromosome 1 is much bigger than chromosome 22.
        cache_jobs: list[tuple[TargetGenome, tuple[tuple[str, int], ...]]] = get_jobs(
            is_vcfs_split_by_chrom=is_vcfs_split_by_chrom,
            target_genomes=target_genomes,
            unique_positions=positions_to_query,
            batch_size=args.batch_size,
        )
        max_workers: int = min(args.workers, len(cache_jobs))
        task_id: TaskID = progress.add_task(
            start=True,
            total=len(cache_jobs),
            description="Caching variants...",
        )

        if max_workers == 1:
            run_jobs_sequentially(progress=progress, jobs=cache_jobs, task_id=task_id)
        else:
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=max_workers,
                initializer=functools.partial(init_worker, log_queue, args.verbose),
            ) as executor:
                try:
                    run_jobs_parallel(
                        progress=progress,
                        jobs=cache_jobs,
                        task_id=task_id,
                        executor=executor,
                        timeout_s=args.timeout_seconds,
                    )
                except Exception:
                    logger.info("Shutting down worker processes")
                    executor.shutdown(wait=False, cancel_futures=True)
                    log_queue.put(None)
                    lp.join()

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
        cache_variants(target=target_genome, positions=positions)
        progress.update(task_id, completed=i, refresh=True)


def run_jobs_parallel(
    progress: Progress,
    jobs: list[tuple[TargetGenome, tuple[tuple[str, int], ...]]],
    task_id: TaskID,
    executor: concurrent.futures.Executor,
    timeout_s: int,
) -> None:
    futures: list[concurrent.futures.Future] = []
    for job in jobs:
        target_genome: TargetGenome
        positions: Sequence[tuple[str, int]]
        target_genome, positions = job
        futures.append(
            executor.submit(
                cache_variants,
                target=target_genome,
                positions=positions,
            )
        )

    try:
        for i, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            _ = future.result(timeout=timeout_s)
            progress.update(task_id, completed=i, refresh=True)
            progress.print(f"Job {i} finished of {len(jobs)} total")
    except Exception as e:
        logger.critical("A worker process has stopped suddenly, bailing out")
        executor.shutdown(cancel_futures=True, wait=True)
        raise e


def cleanup_cache_tmpdir(cache_dir: pathlib.Path) -> None:
    tmp_dir = cache_dir / "tmp"
    try:
        rmtree(tmp_dir)
    except FileNotFoundError:
        pass
