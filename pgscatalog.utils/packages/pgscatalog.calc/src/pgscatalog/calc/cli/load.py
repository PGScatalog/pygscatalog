from __future__ import annotations

import argparse
import atexit
import contextlib
import logging
import logging.handlers
import pathlib
import shutil
import zipfile
from dataclasses import dataclass
from itertools import batched
from shutil import rmtree

import zarr
import zarr.storage
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
from pgscatalog.calc.lib.cache.genotypecache import leaf_group_path
from pgscatalog.calc.lib.config import config
from pgscatalog.calc.lib.constants import canonicalise_chrom

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

        # the samplesheet is a boundary too: `path:chrX` and `path:X` are one leaf
        return cls(path=path, chrom=canonicalise_chrom(chrom))

    def leaf_path(self, sampleset: str) -> str:
        """Where this input lands in the cache hierarchy: the same string LOAD writes

        Parameters
        ----------
        sampleset : str
            The sampleset these inputs are being loaded into.
        """
        return leaf_group_path(
            sampleset=sampleset,
            chromosome=self.chrom,
            filename=pathlib.Path(self.path).name,
        )


def check_leaf_paths(targets: list[AnnotatedGenomePath], sampleset: str) -> None:
    """Check inputs don't collide in the cache hierarchy before launching any work

    Directories are stripped from input paths, so two files with the same basename
    on the same chromosome share a zarr leaf group. LOAD is fanned out per shard,
    so the collision would only appear when the per-shard archives are merged, and
    one leaf would silently overwrite the other (ADR 0003).

    Pure path arithmetic, no file access: the leaf path is the one `TargetGenome`
    writes to, so this guards the string that actually collides.
    """
    seen: dict[str, str] = {}

    for target in targets:
        leaf_path = target.leaf_path(sampleset)
        if (clash := seen.get(leaf_path)) is not None:
            raise ValueError(
                f"Target genomes {clash!r} and {target.path!r} both cache to "
                f"{leaf_path!r}. Directories are stripped from input paths, "
                "so rename one file or stage it with a unique name (stageAs)."
            )
        seen[leaf_path] = target.path


def parse_genome_paths(
    target_genomes: list[str],
    cache_dir: Pathish,
    sampleset: str,
    bgen_sample_file: Pathish | None,
) -> list[TargetGenome]:
    targets: list[AnnotatedGenomePath] = sorted(
        {AnnotatedGenomePath.from_string(x) for x in target_genomes},
        key=lambda x: (x.path, x.chrom or ""),
    )

    check_leaf_paths(targets, sampleset)

    genomes = []

    for genome in targets:
        genomes.append(
            TargetGenome(
                target_path=genome.path,
                target_index_path=None,
                chrom=genome.chrom,
                cache_dir=cache_dir,
                sampleset=sampleset,
                sample_file=bgen_sample_file,
            )
        )

    return genomes


def check_chromosomes(target_genomes: list[TargetGenome]) -> None:
    """Check chromosomes are valid before launching any work"""
    if unsupported := {x.chrom for x in target_genomes} - VALID_CHROMOSOMES:
        names = ", ".join(sorted(str(chrom) for chrom in unsupported))
        raise ValueError(f"{names} not supported yet: 1-22, X, Y and MT can be scored")


def load_cli(args: argparse.Namespace) -> None:
    # override zarr chunk size
    config.ZARR_VARIANT_CHUNK_SIZE = args.batch_size

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

    scorefiles: Scorefiles = Scorefiles(args.score_paths)

    for target_genome in track(
        target_genomes, description="Querying target genome index..."
    ):
        if args.zarr_zip_file is not None:
            group_path = f"genotypes.zarr/{target_genome._zarr_group_path}"
            logger.info(f"--zarr_zip_file set, extracting {group_path=}")
            unzip_zarr(
                zip_path=args.zarr_zip_file,
                cache_path=args.cache_dir,
                group_path=group_path,
            )
            logger.info("Cache set up")

        positions_to_query = scorefiles.get_unique_positions(
            chrom=target_genome.chrom, zarr_group=target_genome.zarr_group
        )

        if len(positions_to_query) == 0:
            print("All variants are already cached. Yay!")
            cleanup_empty_zarr_group(group=target_genome.zarr_group)
            continue

        # extract variants from target genomes in batches of ~5_000
        # build genotype / variant arrays, and then write to zarr
        # (see --batch_size parameter)
        position_batches = list(batched(positions_to_query, args.batch_size))
        for i, positions in enumerate(position_batches, start=1):
            logger.info(
                f"Processing batch {i} of {len(position_batches)} ({args.batch_size=})"
            )
            target_genome.cache_variants(positions)

        # record the declaration on the leaf, not the sampleset: LOAD is fanned out
        # per shard and merging the zips is last-write-wins on an ancestor zarr.json,
        # so a sampleset-level attr would swallow a contradictory declaration instead
        # of letting MATCH detect it (ADR 0004)
        target_genome.zarr_group.attrs["absence_semantics"] = args.absence_semantics

    # prepare to move the zarr array into a read-only single file archive
    logger.info("Packaging zarr zip store")
    zip_path = args.cache_dir / "genotypes.zarr.zip"
    zarr_path = args.cache_dir / "genotypes.zarr"

    # check if any arrays were written to the zarr store
    store = zarr.storage.LocalStore(zarr_path)
    group = zarr.open_group(store=store)

    if not has_any_arrays(group):
        logger.warning(
            "No arrays in zarr store, so not creating the standard zip output. "
            "This is normal if all index queries returned no variants. "
            "If you're not expecting this: do your genotypes have good "
            "variant density? Are they imputed? Are your scorefile positions/build OK?"
        )
        return

    zip_zarr(zip_path=zip_path, zarr_path=zarr_path)

    # clean up the original directory store
    logger.info("Cleaning up zarr directory store")
    shutil.rmtree(zarr_path)

    print("Finished caching :tada: Goodbye!")


def has_any_arrays(group: zarr.Group) -> bool:
    """Does a zarr group contain any arrays? Recurses the hierarchy"""
    for key in group:
        item = group[key]

        if isinstance(item, zarr.Array):
            return True

        if isinstance(item, zarr.Group) and has_any_arrays(item):
            return True

    return False


def cleanup_empty_zarr_group(group: zarr.Group) -> None:
    """
    Delete  a group that doesn't contain a genotype array.

    This can happen if no results were returned after querying an index.

    Empty groups won't have the correct interface for consumers, so delete them
    """
    try:
        _ = group["genotypes"]
    except KeyError:
        root = zarr.open_group(store=group.store)
        logger.info(f"No genotypes array found, deleting {group.path=}")
        del root[group.path]


def _members_for_group(names: list[str], group_path: str) -> list[str]:
    """The archive members that make up one leaf group, and nothing else

    A leaf's own files, plus the `zarr.json` of every ancestor group: without
    those the hierarchy doesn't parse and a downstream stage can't open the leaf
    (ADR 0003).
    """
    parts = group_path.split("/")
    ancestors = {"/".join(parts[:i]) + "/zarr.json" for i in range(1, len(parts))}
    prefix = f"{group_path}/"
    return [n for n in names if n.startswith(prefix) or n in ancestors]


def unzip_zarr(
    zip_path: pathlib.Path, cache_path: pathlib.Path, group_path: str | None = None
) -> None:
    """Extract a zarr zip store to a directory

    This means zarr can work with a DirectoryStore to update an exising cache

    Optionally extract a specific group path from the zip file
    """
    if not zip_path.exists():
        raise FileNotFoundError(
            f"{zip_path=} does not existRemove --zarr_zip_file argument"
        )

    logger.info(f"Extracting zarr zip store to directory {cache_path}")

    with zipfile.ZipFile(zip_path) as archive:
        members = None
        if group_path is not None:
            logger.info(f"Only extracting {group_path=} from zip")
            members = _members_for_group(archive.namelist(), group_path)
        archive.extractall(path=cache_path, members=members)

    logger.info("Unzipping finished")


def zip_zarr(zip_path: pathlib.Path, zarr_path: pathlib.Path) -> None:
    """Archive a zarr directory store

    Stored, not deflated: zarr has already compressed the arrays, so the bytes are
    copied rather than recompressed (ADR 0012).
    """
    logger.info(f"Archiving {zarr_path} to {zip_path}")

    # paths are relative to the store's parent, so the archive holds
    # genotypes.zarr at its top level wherever it's extracted
    root = zarr_path.parent

    with zipfile.ZipFile(zip_path, mode="w", compression=zipfile.ZIP_STORED) as archive:
        for path in sorted(zarr_path.rglob("*")):
            if path.is_file():
                archive.write(path, arcname=path.relative_to(root))

    logger.info("Zipping finished")


def cleanup_cache_tmpdir(cache_dir: pathlib.Path) -> None:
    tmp_dir = cache_dir / "tmp"
    with contextlib.suppress(FileNotFoundError):
        rmtree(tmp_dir)
