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
    MISSING_GENOTYPE_SENTINEL_VALUE,
    MISSING_GENOTYPE_TUPLE,
    NO_ALLELE_SENTINEL_VALUE,
    ZARR_MAX_PLOIDY,
    canonicalise_chrom,
)

from .targetvariants import (
    TargetVariants,
    add_missing_positions_to_lists,
    collapse_mt_calls,
)

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable, Sequence
    from typing import IO

    from numpy import typing as npt

    from pgscatalog.calc.lib.types import Pathish

logger = logging.getLogger(__name__)


def vcf_get_genotypes(
    record_samples: pysam.VariantRecordSamples,
) -> Generator[tuple[int, ...]]:
    """Read ploidy from the length of the GT field, padding to the width of the
    widest call (ADR 0007).

    A hemizygous call keeps its one allele and fills the slot it doesn't have with
    the no-allele sentinel. A wholly missing call is stored as a missing diploid
    call whatever its GT length, so on a haploid chromosome it is imputed and
    counted at ploidy 2 rather than 1. That is the one place this reader departs
    from ADR 0007's encoding table, which gives a hemizygous no-call [255, 254].
    """
    for sample in record_samples.values():
        gt: tuple[int | None, ...] | None = sample.get("GT", None)
        if gt is None or all(allele is None for allele in gt):
            yield MISSING_GENOTYPE_TUPLE
        elif len(gt) > ZARR_MAX_PLOIDY:
            # keeping the first two alleles would turn a call the cache can't hold
            # into a plausible score
            raise ValueError(f"Ploidy {len(gt)} is greater than {ZARR_MAX_PLOIDY=}")
        else:
            alleles = tuple(
                MISSING_GENOTYPE_SENTINEL_VALUE if allele is None else allele
                for allele in gt
            )
            yield alleles + (NO_ALLELE_SENTINEL_VALUE,) * (
                ZARR_MAX_PLOIDY - len(alleles)
            )


def vcf_get_sample_list(target_path: Pathish) -> list[str]:
    with pysam.VariantFile(str(target_path), mode="r") as vcf:
        return list(vcf.header.samples)


def region_sort_key(region: tuple[str, int, int]) -> tuple[int, str, int]:
    """Sort canonical regions by chromosome, then position.

    Canonical names aren't all numeric (ADR 0008), so numbered chromosomes sort
    numerically ahead of X, Y and MT instead of raising on int("X").
    """
    chrom, start, _ = region
    return (int(chrom) if chrom.isdigit() else 99, chrom, start)


def vcf_buffer_variants(
    position_batch: Iterable[tuple[str, int]],
    target_path: Pathish,
    index_path: Pathish,
    cache_dir: Pathish,
    sampleset: str,
    target_chroms: Sequence[str] = (),
) -> TargetVariants:
    """
    Query a target genome with a tabix subprocess, buffering results to a temporary
     VCF file

    The temporary files are read to yield TargetVariants and automatically cleaned up

    The aim of this approach is to minimise RAM usage at the cost of some local storage

    target_chroms are the contigs the VCF header declares. Positions arrive with
    canonical chromosome names, so they're translated back into whatever dialect this
    target uses before bcftools sees them - a chr-prefixed VCF asked for region "1"
    returns nothing at all (ADR 0008).
    """
    regions: list[tuple[str, int, int]] = sorted(
        [(x[0], x[1], x[1]) for x in position_batch], key=region_sort_key
    )
    dialect = {canonicalise_chrom(chrom): chrom for chrom in target_chroms}
    query_regions = [(dialect.get(c, c), start, end) for c, start, end in regions]
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
                regions=query_regions,
                region_file=region_file,
                target_path=target_path,
                index_path=index_path,
                variant_buffer=variant_buffer,
            )
            # read the temporary file and make a TargetVariants class
            return parse_target_variants(
                target_path=variant_buffer.name,
                scoring_file_regions=regions,
                sampleset=sampleset,
            )
    finally:
        pathlib.Path(region_file.name).resolve().unlink(missing_ok=True)
        pathlib.Path(variant_buffer.name).resolve().unlink(missing_ok=True)


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
            # the cache stores canonical names only (ADR 0008)
            chrom = canonicalise_chrom(variant.chrom)
            seen_positions.add((chrom, variant.pos))
            # every call is padded to the same width, so the array stays rectangular
            gts: npt.NDArray[np.uint8] = np.fromiter(
                vcf_get_genotypes(variant.samples),
                dtype=(np.uint8, ZARR_MAX_PLOIDY),
            )

            if chrom == "MT":
                collapse_mt_calls(gts, variant.pos, target_path)

            if chrom == "Y" and np.any(gts == MISSING_GENOTYPE_SENTINEL_VALUE):
                # a male-only chrY file is unambiguous - every call in it is
                # hemizygous - but a `.` is not: VCF writes it for a no-call and for
                # a sample that has no Y alike, and nothing here can tell them apart
                raise ValueError(
                    f"Missing chrY call at Y:{variant.pos} in {target_path}. VCF "
                    "writes `.` for both a no-call and a sample that has no Y, so "
                    "this pipeline can't tell which this is. Reading ploidy needs a "
                    "format that records it, like bgen (ADR 0007)"
                )

            alts = list(variant.alts) if variant.alts is not None else None

            chroms.append(chrom)
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
    index_path: Pathish,
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

    # handle odd index cases https://github.com/samtools/bcftools/issues/984
    # e.g. symlinked to different directories
    target_str = str(pathlib.Path(target_path).resolve())
    index_str = str(pathlib.Path(index_path).resolve())
    target_path_with_index = f"{target_str}##idx##{index_str}"

    pysam.bcftools.view(
        "-R",
        str(region_file.name),
        target_path_with_index,
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
