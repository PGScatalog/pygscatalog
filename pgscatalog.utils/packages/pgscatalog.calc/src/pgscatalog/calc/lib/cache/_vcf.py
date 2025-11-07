from __future__ import annotations

import csv
import logging
import pathlib
import time
from collections.abc import Sequence
from tempfile import NamedTemporaryFile
from typing import TYPE_CHECKING

import numpy as np
import pysam
import pysam.bcftools

from pgscatalog.calc.lib.constants import (
    MISSING_GENOTYPE_TUPLE,
)

from .targetvariants import TargetVariants, add_missing_positions_to_lists

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable, Sequence
    from typing import IO

    from numpy import typing as npt

    from pgscatalog.calc.lib.types import Pathish

logger = logging.getLogger(__name__)


def vcf_get_genotypes(
    record_samples: pysam.VariantRecordSamples,
) -> Generator[tuple[int, int]]:
    for sample in record_samples.values():
        gt: tuple[int, int] | None = sample.get("GT", None)
        if gt == (None,) or gt is None:
            yield MISSING_GENOTYPE_TUPLE
        else:
            yield gt


def vcf_get_sample_list(target_path: Pathish) -> list[str]:
    with pysam.VariantFile(str(target_path), mode="r") as vcf:
        return list(vcf.header.samples)


def vcf_buffer_variants(
    position_batch: Iterable[tuple[str, int]],
    target_path: Pathish,
    cache_dir: Pathish,
    sampleset: str,
) -> TargetVariants:
    """
    Query a target genome with a tabix subprocess, buffering results to a temporary
     VCF file

    The temporary files are read to yield TargetVariants and automatically cleaned up

    The aim of this approach is to minimise RAM usage at the cost of some local storage
    """
    regions: list[tuple[str, int, int]] = sorted(
        [(x[0], x[1], x[1]) for x in position_batch],
        key=lambda t: (int(t[0]), t[1]),
    )
    tmp_dir = pathlib.Path(cache_dir).resolve() / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # temporary files will get deleted when context manager exits
    try:
        with (
            NamedTemporaryFile(
                dir=tmp_dir, mode="wt", delete_on_close=False
            ) as region_file,
            NamedTemporaryFile(
                dir=tmp_dir, delete_on_close=False, suffix=".bcf"
            ) as variant_buffer,
        ):
            # query the indexed VCF, creating a temporary file
            run_bcftools_view(
                regions=regions,
                region_file=region_file,
                target_path=target_path,
                variant_buffer=variant_buffer,
            )
            # read the temporary file and make a TargetVariants class
            return parse_target_variants(
                target_path=variant_buffer.name,
                scoring_file_regions=regions,
                sampleset=sampleset,
            )
    finally:
        pathlib.Path(region_file.name).unlink()
        pathlib.Path(variant_buffer.name).unlink()


def parse_target_variants(
    target_path: Pathish,
    scoring_file_regions: Sequence[tuple[str, int, int]],
    sampleset: str,
) -> TargetVariants:
    start = time.perf_counter()
    seen_positions = set()
    chroms: list[str] = []
    positions: list[int] = []
    ref_alleles: list[str | None] = []
    alt_alleles: list[list[str] | None] = []
    hard_calls: list[npt.NDArray[np.uint8]] = []

    with pysam.VariantFile(str(target_path), mode="r") as records:
        samples = list(records.header.samples)
        n_samples = len(samples)

        i = 0
        for i, variant in enumerate(records.fetch()):
            seen_positions.add((variant.chrom, variant.pos))
            try:
                gt_tuples: Iterable[tuple[int, int]] = list(
                    vcf_get_genotypes(variant.samples)
                )
                gts: npt.NDArray[np.uint8] = np.fromiter(gt_tuples, dtype=(np.uint8, 2))
            except ValueError as e:
                # happens when some samples are haploid
                raise NotImplementedError(
                    "Sex chromosome genotypes not supported"
                ) from e

            alts = list(variant.alts) if variant.alts is not None else None

            chroms.append(variant.chrom)
            positions.append(variant.pos)
            ref_alleles.append(variant.ref)
            alt_alleles.append(alts)
            hard_calls.append(gts)

        end = time.perf_counter()
        logger.info(
            f"Processed {i / (end - start):.2f} variants per second in local buffer"
        )

        # add any unseen variants to the lists (mutating in place!)
        add_missing_positions_to_lists(
            chroms=chroms,
            positions=positions,
            ref_alleles=ref_alleles,
            alt_alleles=alt_alleles,
            hard_calls=hard_calls,
            scoring_file_regions=[
                (chrom, pos) for chrom, pos, _ in scoring_file_regions
            ],
            seen_positions=seen_positions,
            n_samples=n_samples,
        )

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


def run_bcftools_view(
    regions: list[tuple[str, int, int]],
    target_path: Pathish,
    variant_buffer: IO[bytes],
    region_file: IO[str],
) -> None:
    logger.info(f"Buffering {len(regions)=} variants to temporary file")
    csv.writer(region_file, delimiter="\t").writerows(regions)
    region_file.close()  # closed but used by bcftools view

    logger.info(f"Running bcftools view: {target_path} saved to {variant_buffer.name}")
    start = time.perf_counter()
    # variant buffer is a temporary uncompressed bcf file
    # --exclude-types: excludes symbolic alleles (see man bcftools, basically not
    #   snp / mnp / indel)
    # --exclude-uncalled: only include sites with at least one non-missing genotype
    # --write-index: suppress annoying htslib warnings when working with the file
    pysam.bcftools.view(
        "-R",
        str(region_file.name),
        str(pathlib.Path(target_path).resolve()),
        "-o",
        str(variant_buffer.name),
        "--output-type",
        "u",
        "--exclude-types",
        "other",
        "--exclude-uncalled",
        "--write-index",
        # don't ingest multiallelics in RC1
        "--min-alleles",
        "2",
        "--max-alleles",
        "2",
        catch_stdout=False,
    )
    end = time.perf_counter()
    logger.info(
        f"bcftools view finished in {end - start:.2f} seconds: "
        f"{len(regions) / (end - start):.2f} positions per second"
    )
