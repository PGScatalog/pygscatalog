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

from bgen import BgenReader

from .targetvariant import TargetVariant

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Iterable, Sequence

    from .types import Pathish


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


def parse_target_variants(
    target_path: Pathish,
    scoring_file_regions: Sequence[tuple[str, int]],
    skip_multiallelic: bool = True,
    ref_first: bool = True,
) -> Generator[TargetVariant, None, None]:
    seen_positions = set()
    start = time.perf_counter()
    with BgenReader(str(target_path), delay_parsing=True) as bgenfile:
        logger.info(f"{target_path=} opened")
        i = 0
        for variant in bgenfile:
            i += 1
            # ignore multiallelics for now
            if len(variant.alleles) > 2 and skip_multiallelic:
                logger.info(f"Skipping multiallelic variant {variant=}")
                continue
            elif len(variant.alleles) > 2 and not skip_multiallelic:
                raise NotImplementedError("Can't handle multiallelic variants")

            # cast list to tuple for @cache'd is_snp function (input must be hashable)
            if not is_snp(alleles=tuple(variant.alleles)):
                logger.info(f"Skipping non-SNP {variant.alleles=}")
                continue

            # use PGS Catalog style chromosomes (no padding or prefix) to build cache
            fetched_chrom = str(variant.chrom).lstrip("chr").lstrip("0")
            fetched_pos = variant.pos
            seen_positions.add((fetched_chrom, fetched_pos))

            ref: str = variant.alleles[0] if ref_first else variant.alleles[1]
            alt: str = variant.alleles[1] if ref_first else variant.alleles[0]

            # if variant.is_phased:
            #     gts = phased_probabilities_to_hard_calls(variant.probabilities)
            # else:
            #     gts = unphased_probabilities_to_hard_calls(variant.probabilities)

            yield TargetVariant(
                chr_name=fetched_chrom,
                chr_pos=fetched_pos,
                ref=ref,
                alts=(alt,),
                probs=variant.probabilities,
                is_phased=variant.is_phased,
            )

    end = time.perf_counter()
    logger.info(
        f"Processed {i / (end - start):.2f} variants per second in local buffer"
    )

    missing_positions = set(scoring_file_regions) - seen_positions
    for chr_name, chr_pos in missing_positions:
        # important to prevent future cache misses
        logger.debug(f"{chr_name=} {chr_pos=} not found in in {target_path=}")
        yield TargetVariant(
            chr_name=chr_name, chr_pos=chr_pos, ref=None, alts=None, gts=None
        )


def bgen_buffer_variants(
    cache_dir: Pathish,
    position_batch: Sequence[tuple[str, int]],
    target_path: Pathish,
    idx_path: Pathish,
    target_chroms: list[str],
    skip_multiallelic: bool = True,
    ref_first: bool = True,
) -> Iterable[TargetVariant]:
    if shutil.which("bgenix") is None:
        raise FileNotFoundError("Error: bgenix program not found in PATH")

    tmp_dir = pathlib.Path(cache_dir).resolve() / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    buffered_bgen_file = cast(
        BinaryIO,
        NamedTemporaryFile(
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
        # close the temporary bgen file and get ready to read it
        buffered_bgen_file.close()

        # parse bgen variants as target variants
        yield from parse_target_variants(
            target_path=buffered_bgen_path,
            scoring_file_regions=position_batch,
            skip_multiallelic=skip_multiallelic,
            ref_first=ref_first,
        )
    finally:
        buffered_bgen_path.unlink()
        with sqlite3.connect(idx_path) as conn:
            if temp_table_name is not None:
                conn.execute(f"DROP TABLE {temp_table_name}")
            if temp_view_name is not None:
                conn.execute(f"DROP VIEW {temp_view_name}")


def create_bgen_index_query(
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
