from __future__ import annotations

import csv
import logging
import pathlib
import time
from tempfile import NamedTemporaryFile
from typing import TYPE_CHECKING

import numpy as np
import pysam
import pysam.bcftools

from .constants import MISSING_GENOTYPE_TUPLE

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable

    from numpy import typing as npt

    from .types import Pathish

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
):
    """
    Query a target genome with a tabix subprocess, buffering results to a temporary
     VCF file

    The temporary files are read to yield TargetVariants and automatically cleaned up

    The aim of this approach is to minimise RAM usage at the cost of some local storage
    """
    # tabix regions file format is CHROM \t START \t END (1-based)
    try:
        regions: list[tuple[str, int, int]] = sorted(
            [(x[0], x[1], x[1]) for x in position_batch],
            key=lambda t: (int(t[0]), t[1]),
        )
    except ValueError as e:
        raise NotImplementedError("Only chromosomes 1 -22 supported") from e

    tmp_dir = pathlib.Path(cache_dir).resolve() / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # temporary files will get deleted when context manager exits
    with (
        NamedTemporaryFile(
            dir=tmp_dir, mode="wt", delete_on_close=False
        ) as region_file,
        NamedTemporaryFile(
            dir=tmp_dir, delete_on_close=False, suffix=".bcf"
        ) as variant_buffer,
    ):
        logger.info(f"Buffering {len(regions)=} variants to temporary file")
        csv.writer(region_file, delimiter="\t").writerows(regions)
        region_file.close()  # closed but used by bcftools view

        logger.info(
            f"Running bcftools view: {target_path} saved to {variant_buffer.name}"
        )
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

        seen_positions = set()
        start = time.perf_counter()

        with pysam.VariantFile(variant_buffer.name, mode="r") as records:
            i = 0
            for variant in records.fetch():
                seen_positions.add((variant.chrom, variant.pos))
                try:
                    gt_tuples: Iterable[tuple[int, int]] = list(
                        vcf_get_genotypes(variant.samples)
                    )
                    gts: npt.NDArray[np.uint8] = np.fromiter(
                        gt_tuples, dtype=(np.uint8, 2)
                    )
                except ValueError as e:
                    # happens when some samples are haploid
                    raise NotImplementedError(
                        "Sex chromosome genotypes not supported"
                    ) from e
                else:
                    yield TargetVariant(
                        chr_name=variant.chrom,
                        chr_pos=variant.pos,
                        ref=variant.ref,
                        alts=variant.alts,
                        gts=gts,
                    )
                    i += 1

        end = time.perf_counter()
        logger.info(
            f"Processed {i / (end - start):.2f} variants per second in local buffer"
        )

        missing_positions = set(position_batch) - seen_positions
        for chr_name, chr_pos in missing_positions:
            # important to prevent future cache misses
            logger.debug(f"{chr_name=} {chr_pos=} not found in in {target_path=}")
            yield TargetVariant(
                chr_name=chr_name, chr_pos=chr_pos, ref=None, alts=None, gts=None
            )
