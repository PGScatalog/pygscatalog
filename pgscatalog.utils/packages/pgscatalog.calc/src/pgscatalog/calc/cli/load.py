from __future__ import annotations

import argparse
import atexit
import contextlib
import logging
import logging.handlers
import shutil
import subprocess
from dataclasses import dataclass
from shutil import rmtree
from typing import TYPE_CHECKING

from rich import print
from rich.logging import RichHandler
from rich.progress import track

from pgscatalog.calc import (
    VALID_CHROMOSOMES,
    GenomeFileType,
    Pathish,
    Scorefiles,
    TargetGenome,
)

if TYPE_CHECKING:
    import pathlib

logger = logging.getLogger(__name__)

if shutil.which("7z") is None:
    raise ValueError(
        "7z application not detected on your system."
        "Please install 7z or use bioconda to install this package"
    )

if shutil.which("bgenix") is None:
    raise ValueError(
        "bgenix application not detected on your system."
        "Please install bgenix or use bioconda to install this package."
    )


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
    """Check chromosomes are valid before launching any work"""
    if not {x.chrom for x in target_genomes}.issubset(VALID_CHROMOSOMES):
        raise NotImplementedError("Chromosomes 1-22 supported but this will be fixed")

    is_chrom_split = all(x.chrom is not None for x in target_genomes)

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

    print(":test_tube:", "pgsc_calc load", ":dna:", "starting...")

    # set up input data
    target_genomes: list[TargetGenome] = parse_genome_paths(
        target_genomes=args.target_genomes,
        cache_dir=args.cache_dir,
        sampleset=args.sampleset,
        bgen_sample_file=args.bgen_sample_file,
    )

    # check chromosome input
    check_chromosomes(target_genomes=target_genomes)

    # has an existing cache been provided?
    if args.zarr_zip_file is not None:
        logger.info("--zarr_zip_file set, updating existing cache")
        unzip_zarr(zip_path=args.zarr_zip_file, cache_path=args.cache_dir)
    else:
        logger.info("--zarr_zip_file not set, creating new cache")

    scorefiles: Scorefiles = Scorefiles(args.score_paths)

    for target_genome in track(
        target_genomes, description="Querying target genome index..."
    ):
        positions_to_query = scorefiles.get_unique_positions(
            chrom=target_genome.chrom, zarr_group=target_genome.zarr_group
        )

        if len(positions_to_query) == 0:
            print("All variants are already cached. Yay!")
            continue
        target_genome.cache_variants(positions_to_query)

    # move the zarr array into a read-only single file archive
    logger.info("Packaging zarr zip store")
    zip_path = args.cache_dir / "genotypes.zarr.zip"
    zarr_path = args.cache_dir / "genotypes.zarr"
    zip_zarr(zip_path=zip_path, zarr_path=zarr_path)

    # clean up the original directory store
    logger.info("Cleaning up zarr directory store")
    shutil.rmtree(zarr_path)

    print("Finished caching :tada: Goodbye!")


def unzip_zarr(zip_path: pathlib.Path, cache_path: pathlib.Path) -> None:
    """Run a 7z subprocess to unzip a zarr store

    This means zarr can work with a DirectoryStore to update an exising cache
    """
    if not zip_path.exists():
        raise FileNotFoundError(
            f"{zip_path=} does not existRemove --zarr_zip_file argument"
        )

    if (zarr_archive := cache_path / "genotypes.zarr").exists():
        raise FileExistsError(
            f"{zarr_archive=} already existsChoose a different --cache_dir"
        )

    logger.info(f"Extracting zarr zip store to directory {cache_path}")
    cmd = ["7z", "x", str(zip_path), f"-o{str(cache_path)}"]
    logger.info(f"Running subprocess {cmd=}")
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)

    if result.stdout:
        logger.info(result.stdout)

    if result.stderr:
        logger.error(result.stderr)

    logger.info("Unzipping finished")


def zip_zarr(zip_path: pathlib.Path, zarr_path: pathlib.Path) -> None:
    """Run a 7z subprocess to create a zip archive

    zipfile in standard library and other python libraries are pure Python. 7z is
    written in C++ so shouldn't have problems scaling.
    """
    cmd = [
        "7z",
        "a",
        "-tzip",  # make a zip archive
        "-mx0",  # don't compress (zarr already compresses arrays)
        str(zip_path),
        str(zarr_path),
    ]

    logger.info(f"Running subprocess {cmd=}")
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)

    if result.stdout:
        logger.info(result.stdout)

    if result.stderr:
        logger.error(result.stderr)

    logger.info("Zipping finished")


def cleanup_cache_tmpdir(cache_dir: pathlib.Path) -> None:
    tmp_dir = cache_dir / "tmp"
    with contextlib.suppress(FileNotFoundError):
        rmtree(tmp_dir)
