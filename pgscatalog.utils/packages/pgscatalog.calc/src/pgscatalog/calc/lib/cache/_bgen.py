"""
Tools for extracting and converting variants from indexed BGEN files.

This module provides functions to read, filter, and normalise variant data stored in
BGEN format, typically used for large-scale imputed genotype datasets. The internal
TargetVariants class is returned to provide a consistent interface.

Key functionality includes:

    Hard call conversion: Convert phased and unphased genotype probabilities into
    hard-called genotypes (see phased_to_hardcall and unphased_to_hardcall).

    Sample extraction: Retrieve the list of sample identifiers from a BGEN file.

    Chromosome normalisation: BGEN files can encode chromosomes with a chr prefix or
    0-padding. Queries are translated from PGS Catalog style to match the target genome
    format.

    bgenix subprocess: bgenix is used to create new temporary bgen files. The new bgen files
    get iterated over. This is motivated to promote batch caching, and allow checkpointing
    in case jobs fail. To avoid bgenix query limitations, the index file (a SQLite database)
    is copied and modified.

    Variant parsing: Transform BGEN variants in a temporary bgen file into a validated
    TargetVariants object suitable for downstream data pipelines.

External dependencies:

    bgen: for reading BGEN files directly from Python.
    bgenix: an external command-line utility for querying indexed BGEN files.
    numpy: for genotype probability and hard-call array handling

Notes:

    Only SNP variants are supported, non-SNPs are ignored
    Missing genotypes are represented by a sentinel value defined in constants.py
    Multiallelic variants are not yet supported but internal data structures expect
    a list of alternate alleles
"""

from __future__ import annotations

import logging
import pathlib
import shutil
import sqlite3
import subprocess
import time
from functools import cache
from tempfile import NamedTemporaryFile
from typing import TYPE_CHECKING, BinaryIO, cast

import numpy as np
import numpy.typing as npt
from bgen import BgenReader  # type: ignore[import-untyped]

from ..constants import MISSING_GENOTYPE_SENTINEL_VALUE
from .targetvariants import TargetVariants

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Sequence

    from ..types import Pathish


GT_LOOKUP_TABLE = np.array(
    [
        [0, 0],
        [0, 1],
        [1, 1],
        [MISSING_GENOTYPE_SENTINEL_VALUE, MISSING_GENOTYPE_SENTINEL_VALUE],
    ],
    dtype=np.uint8,
)


def unphased_to_hardcall(probs: npt.NDArray[np.float64]) -> npt.NDArray[np.uint8]:
    """
    Convert unphased probabilities to hard call genotypes, like in a VCF file

    Args:
        probs: a numpy array of shape (n_samples, 3)

    Returns:
        A numpy array of shape (n_samples, 2), e.g. sample A -> 0/0

    Notes:
         Columns in probs array correspond to homozygous first allele (AA),
         heterozygous (Aa), and homozygous second allele (aa).
    """
    return GT_LOOKUP_TABLE[np.argmax(probs, axis=1).astype(np.uint8)]


def phased_to_hardcall(probs: npt.NDArray[np.float64]) -> npt.NDArray[np.uint8]:
    """
    Convert phased probabilities to hard call genotypes, like in a VCF file

    Args:
        probs: a numpy array of shape (n_samples, 4), which represent a haplotype pair

    Returns:
        A numpy array of shape (n_samples, 2), e.g. sample A -> 0/0

    Notes:
       Probability array: first two for haplotype 1 (hap1-allele1, hap1-allele2),
        last two for haplotype 2 (hap2-allele1, hap2-allele2).
    """
    hap1 = np.argmax(probs[:, :2], axis=1).astype(np.uint8)
    hap2 = np.argmax(probs[:, 2:], axis=1).astype(np.uint8)
    return np.column_stack((hap1, hap2))


def bgen_get_sample_list(target_path: Pathish, sample_path: Pathish) -> list[str]:
    with BgenReader(str(target_path), str(sample_path), delay_parsing=True) as bfile:
        return list(bfile.samples)


def prepare_region_query(
    regions: Iterable[tuple[str, int]], target_chroms: list[str]
) -> list[tuple[str, int, int]]:
    logger.info(f"Preparing regions for bgenix query, {target_chroms=}")

    # set up transform function for chromosome
    transform: Callable
    if any(x.startswith("0") for x in target_chroms):
        # are they padded? e.g. UK Biobank
        logger.info("Padded chromosomes detected")
        transform = add_chrom_padding
    elif any(x.startswith("chr") for x in target_chroms):
        logger.info("Prefixed chromosones detected")
        # are they prefixed? e.g. some imputation servers
        transform = add_chrom_prefix
    else:
        logger.info("PGS Catalog format chromosomes detected")
        transform = lambda x: x  # noqa: E731

    return [(transform(chrom), pos, pos) for chrom, pos in regions]


def run_bgenix_subprocess(
    *,
    target_path: Pathish,
    temp_bgen_file: BinaryIO,
    temp_idx_path: pathlib.Path,
    table_name: str,
) -> None:
    cmd = [
        "bgenix",
        "-g",
        str(target_path),
        "-i",
        str(temp_idx_path),
        "-table",
        table_name,
    ]
    logger.info(f"Running bgenix in a subprocess {cmd=}")
    start = time.perf_counter()
    bgenix_output = subprocess.run(
        cmd,
        stdout=temp_bgen_file,
        stderr=subprocess.DEVNULL,
    )
    end = time.perf_counter()
    logger.info(f"bgenix finished in {end - start:.2f} seconds.")

    if bgenix_output.returncode != 0:
        raise ValueError("bgenix exited with non-zero return code")


@cache
def is_snp(alleles: tuple[str]) -> bool:
    # squish a list of alleles into a single string and check they're all ACTG
    return set("".join(alleles)).issubset({"A", "C", "G", "T"})


@cache
def normalise_chrom(chrom: str | int) -> str:
    """
    Some bgen files encode chromosomes with a leading 0 or chr prefix
    """
    s = str(chrom)
    if s.startswith("chr"):
        s = s[3:]
    return s.lstrip("0")


def parse_target_variants(
    *,
    target_path: Pathish,
    buffered_path: Pathish,
    bgen_sample_path: Pathish,
    scoring_file_regions: Sequence[tuple[str, int]],
    sampleset: str,
    skip_multiallelic: bool = True,
    ref_first: bool = True,
) -> TargetVariants:
    seen_positions = set()
    start = time.perf_counter()

    chroms: list[str] = []
    positions: list[int] = []
    ref_alleles: list[str | None] = []
    alt_alleles: list[list[str] | None] = []
    hard_calls: list[npt.NDArray[np.uint8]] = []

    with BgenReader(
        str(buffered_path), delay_parsing=True, sample_path=str(bgen_sample_path)
    ) as bgenfile:
        samples: list[str] = bgenfile.samples
        n_samples = len(samples)
        logger.info(f"{buffered_path=} opened, {n_samples} present")

        i = 0  # initialise enumerate in case there's no variants in the file
        for i, variant in enumerate(bgenfile):
            # cast list to tuple for @cache'd is_snp function (input must be hashable)
            if not is_snp(alleles=tuple(variant.alleles)):
                logger.info(f"Skipping non-SNP {variant.alleles=}")
                continue

            is_multiallelic = len(variant.alleles) > 2
            if is_multiallelic and skip_multiallelic:
                logger.info(f"Skipping multiallelic variant {variant=}")
                continue
            if is_multiallelic and not skip_multiallelic:
                raise NotImplementedError("Can't handle multiallelic variants")

            # use PGS Catalog style chromosomes (no padding or prefix) to build cache
            fetched_chrom = normalise_chrom(variant.chrom)
            fetched_pos = variant.pos
            ref: str = variant.alleles[0] if ref_first else variant.alleles[1]
            # store as a list for future multiallelic support
            alt: list[str] = [variant.alleles[1]] if ref_first else [variant.alleles[0]]
            seen_positions.add((fetched_chrom, fetched_pos))

            # convert probabilities to uint8 genotypes ASAP
            if variant.is_phased:
                gts = phased_to_hardcall(variant.probabilities)
            else:
                gts = unphased_to_hardcall(variant.probabilities)

            # update all the lists
            chroms.append(fetched_chrom)
            positions.append(fetched_pos)
            ref_alleles.append(ref)
            alt_alleles.append(alt)
            hard_calls.append(gts)

    end = time.perf_counter()
    logger.info(
        f"Processed {i / (end - start):.2f} variants per second in local buffer"
    )

    # represent missing genotypes with a sentinel value of the same shape
    missing_gts = np.full(
        shape=(n_samples, 2), fill_value=MISSING_GENOTYPE_SENTINEL_VALUE, dtype=np.uint8
    )
    missing_positions = set(scoring_file_regions) - seen_positions
    for chr_name, chr_pos in missing_positions:
        # important to prevent future cache misses
        chroms.append(chr_name)
        positions.append(chr_pos)
        ref_alleles.append(None)
        alt_alleles.append(None)
        hard_calls.append(missing_gts)

    return TargetVariants(
        chr_name=chroms,
        pos=positions,
        refs=ref_alleles,
        alts=alt_alleles,
        gts=hard_calls,
        samples=samples,
        target_path=target_path,
        sampleset=sampleset,
    )


def bgen_buffer_variants(
    *,
    cache_dir: Pathish,
    position_batch: Sequence[tuple[str, int]],
    target_path: Pathish,
    sample_path: Pathish,
    sampleset: str,
    idx_path: Pathish,
    target_chroms: list[str],
    skip_multiallelic: bool = True,
    ref_first: bool = True,
) -> TargetVariants:
    if shutil.which("bgenix") is None:
        raise FileNotFoundError("Error: bgenix program not found in PATH")

    tmp_dir = pathlib.Path(cache_dir).resolve() / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    buffered_bgen_file = cast(
        "BinaryIO",
        NamedTemporaryFile(  # noqa: SIM115
            dir=tmp_dir, delete_on_close=False, suffix=".bgen", delete=False
        ),
    )
    buffered_bgen_path = pathlib.Path(buffered_bgen_file.name).resolve()

    temp_table_name = None
    temp_view_name = None

    try:
        # update PGS Catalog chromosomes to match odd bgen encoding
        regions = prepare_region_query(
            regions=position_batch, target_chroms=target_chroms
        )

        # update the index to include a potentially big query range
        temp_table_name, temp_view_name = create_bgen_index_query(
            idx_path=idx_path,
            temp_target_path=buffered_bgen_path,
            regions=regions,
        )

        # extract subset of variants with bgenix using the newly created table
        run_bgenix_subprocess(
            target_path=target_path,
            temp_bgen_file=buffered_bgen_file,
            temp_idx_path=pathlib.Path(idx_path),
            table_name=temp_view_name,
        )

        # parse bgen variants as target variants
        return parse_target_variants(
            buffered_path=buffered_bgen_path,
            scoring_file_regions=position_batch,
            skip_multiallelic=skip_multiallelic,
            ref_first=ref_first,
            bgen_sample_path=sample_path,
            target_path=target_path,
            sampleset=sampleset,
        )
    finally:
        buffered_bgen_file.close()
        buffered_bgen_path.unlink()
        with sqlite3.connect(idx_path) as conn:
            if temp_table_name is not None:
                conn.execute(f"DROP TABLE {temp_table_name}")
            if temp_view_name is not None:
                conn.execute(f"DROP VIEW {temp_view_name}")


def create_bgen_index_query(
    *,
    idx_path: Pathish,
    temp_target_path: Pathish,
    regions: Sequence[tuple[str, int, int]],
) -> tuple[str, str]:
    """
    Create a new bgen index from a copy of the existing index and add a new table
    containing regions of interest. A bgen index is just a sqlite3 database (yay!).

    bgenix can't process more than 9998 position ranges because of sqlite3 limitations:

    https://enkre.net/cgi-bin/code/bgen/tktview/ac405ee793

    Workaround is to create a new table, populate it, and set the table name at query
    time. A copy is made to prevent conflicts when processing in parallel.
    """
    query_table_name = f"{pathlib.Path(temp_target_path).stem}_Range"
    query_view_name = f"{pathlib.Path(temp_target_path).stem}_RangeView"

    with sqlite3.connect(str(idx_path)) as conn:
        logger.info(f"Creating new table {query_table_name} in {idx_path}")
        # enable write ahead logging to mitigate multiple worker processes
        conn.execute(f"""
        CREATE TABLE {query_table_name}
        (
           chromosome TEXT NOT NULL,
           start INT NOT NULL,
           end INT NOT NULL
        );
        """)
        logger.info(f"Inserting {len(regions)} ranges into QueryRange")
        conn.executemany(
            f"""
        INSERT INTO {query_table_name} VALUES(?, ?, ?);
        """,
            regions,
        )
        conn.commit()

        logger.info(f"Creating {query_view_name} for bgenix CLI (-table arg)")
        conn.execute(f"""
            CREATE VIEW {query_view_name} AS SELECT V.*
            FROM {query_table_name} M
            INNER JOIN Variant V
            WHERE V.chromosome == M.chromosome
            AND V.position BETWEEN M.start AND M.end;
        """)

    return query_table_name, query_view_name


@cache
def add_chrom_padding(chrom: str | int) -> str:
    """
    Return a chromosome value zero-padded to two digits if numeric.

    Parameters
    ----------
    chrom : str or int
        Chromosome identifier from a PGS Catalog scoring file.

    Returns
    -------
    str
        Zero-padded chromosome number (e.g., "01", "02") or the original
        non-numeric chromosome string (e.g., "X", "XY", "MT").
    """
    try:
        return f"{int(chrom):02}"
    except ValueError:
        return str(chrom)  # e.g, X, XY, MT


@cache
def add_chrom_prefix(chrom: str | int) -> str:
    """
    Return a chromosome value with the 'chr' prefix.

    Parameters
    ----------
    chrom : str or int
        Chromosome identifier from a PGS Catalog scoring file.

    Returns
    -------
    str
        Chromosome value with 'chr' prefix (e.g., "chr1", "chrX").
    """
    return f"chr{chrom}"
