"""This module contains a data processing pipeline to format ScoreVariants in a
standard way. Each step in the data processing pipeline is a generator that operates
on a list of ScoreVariants and yields updated ScoreVariants. This makes it easy to
plug in extra steps where needed, and lazily works on millions of objects."""

import logging
import pathlib

import pyliftover

from pgscatalog.core.lib.genomebuild import GenomeBuild
from pgscatalog.core.lib.models import Allele
from pgscatalog.core.lib.pgsexceptions import LiftoverError, EffectTypeError
from pgscatalog.core.lib.effecttype import EffectType

logger = logging.getLogger(__name__)


def normalise(
    scoring_file, drop_missing=False, liftover=False, chain_dir=None, target_build=None
):
    """Order of steps is important:

    1. liftover non-harmonised data (quite rare), failed lifts get None'd
    2. remap harmonised data, failed harmonisations get None'd
    3. log and optionally drop bad variants
    """
    logger.info(
        f"Normalise parameters: {drop_missing=}, {liftover=}, {chain_dir=}, {target_build=}"
    )

    if liftover:
        variants = lift(
            scoring_file=scoring_file,
            harmonised=scoring_file.is_harmonised,
            current_build=scoring_file.genome_build,
            target_build=target_build,
            chain_dir=chain_dir,
        )
    else:
        variants = scoring_file.variants

    variants = remap_harmonised(variants, scoring_file.is_harmonised, target_build)

    if drop_missing:
        variants = drop_hla(variants)

    variants = check_effect_type(variants)
    variants = assign_other_allele(variants)
    variants = check_effect_allele(variants, drop_missing)
    variants = detect_complex(variants)

    if scoring_file.is_wide:
        # wide data must be sorted because check_duplicates requires sorted input
        variants = (x for x in sorted(variants, key=lambda x: x["accession"]))

    variants = check_duplicates(variants)

    return variants


def check_effect_type(variants):
    """Check for non-additive variants and complain if found"""
    for variant in variants:
        if variant.effect_type == EffectType.NONADDITIVE:
            raise EffectTypeError(f"{variant.variant_id=} has unsupported effect type")
        yield variant


def check_duplicates(variants):
    """Check if a scoring file contains multiple variants with the same ID
    ID = chr:pos:effect_allele:other_allele
    """
    seen_ids = {}
    current_accession = None
    n_duplicates = 0
    n_variants = 0
    for variant in variants:
        accession = variant.accession

        if accession != current_accession:
            seen_ids = {}
            current_accession = accession

        if variant.variant_id in seen_ids:
            variant.is_duplicated = True
            n_duplicates += 1

        seen_ids[variant.variant_id] = True

        yield variant
        n_variants += 1

    if n_duplicates > 0:
        logger.warning(
            f"{n_duplicates} of {n_variants} variants are duplicated in: {current_accession}"
        )


def drop_hla(variants):
    """Drop HLA alleles from a list of ScoreVariants

    >>> from pgscatalog.core.lib.models import ScoreVariant
    >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0, "chr_name": "1", "chr_position": 1})
    >>> list(drop_hla([variant])) # doctest: +ELLIPSIS
    [ScoreVariant(..., effect_allele=Allele(allele='A', is_snp=True), ...

    >>> variant = ScoreVariant(**{"effect_allele": "P", "effect_weight": 5, "accession": "test", "row_nr": 0, "chr_name": "1", "chr_position": 1})
    >>> list(drop_hla([variant]))
    []
    """
    n_dropped = 0
    p = Allele(allele="P")
    n = Allele(allele="N")

    for variant in variants:
        match variant:
            case _ if variant.effect_allele in (p, n):
                n_dropped += 1
                continue
            case _:
                yield variant

    logger.warning(f"{n_dropped} HLA alleles detected and dropped")


def assign_other_allele(variants):
    """Check if there's more than one possible other allele, remove if true

    >>> from pgscatalog.core.lib.models import ScoreVariant
    >>> variant = ScoreVariant(**{"chr_position": 1, "rsID": None, "chr_name": "1", "effect_allele": "A", "effect_weight": 5, "other_allele": "A", "row_nr": 0, "accession": "test"})
    >>> list(assign_other_allele([variant]))[0] # doctest: +ELLIPSIS
    ScoreVariant(..., effect_allele=Allele(allele='A', is_snp=True), other_allele=Allele(allele='A', is_snp=True), ...)
    >>> variant = ScoreVariant(**{"chr_position": 1, "rsID": None, "chr_name": "1", "effect_allele": "A", "effect_weight": 5, "other_allele": "A/C", "row_nr": 0, "accession": "test"})
    >>> list(assign_other_allele([variant]))[0] # doctest: +ELLIPSIS
    ScoreVariant(..., effect_allele=Allele(allele='A', is_snp=True), other_allele=None, ...)
    """
    n_dropped = 0
    for variant in variants:
        if getattr(variant.other_allele, "has_multiple_alleles", False):
            n_dropped += 1
            variant.other_allele = None
        yield variant

    if n_dropped > 0:
        logger.warning(f"Multiple other_alleles detected in {n_dropped} variants")
        logger.warning("Other allele for these variants is set to missing")


def remap_harmonised(variants, harmonised, target_build):
    """
    Overwrite key attributes with harmonised data, if available.

    In this case chr_name, chr_position, and other allele are missing.
    Perhaps authors submitted rsID and effect allele originally:

    >>> from pgscatalog.core.lib.models import ScoreVariant
    >>> variant = ScoreVariant(**{"chr_position": 1, "rsID": None, "chr_name": "2", "effect_allele": "A", "effect_weight": 5, "accession": "test", "hm_chr": "1", "hm_pos": 100, "hm_rsID": "testrsid", "hm_inferOtherAllele": "A", "row_nr": 0})
    >>> variant
    ScoreVariant(..., effect_allele=Allele(allele='A', is_snp=True), other_allele=None, ...

    >>> list(remap_harmonised([variant], harmonised=True, target_build=GenomeBuild.GRCh38)) # doctest: +ELLIPSIS
    [ScoreVariant(..., effect_allele=Allele(allele='A', is_snp=True), other_allele=Allele(allele='A', is_snp=True), ...
    """
    if harmonised:
        for variant in variants:
            # using the harmonised field in the header to make sure we don't accidentally overwrite
            # positions with empty data (e.g. in an unharmonised file)
            # if harmonisation has failed we _always_ want to use that information
            variant.chr_name = variant.hm_chr
            variant.chr_position = variant.hm_pos
            if variant.other_allele is None:
                variant.other_allele = variant.hm_inferOtherAllele
            # update the accession to reflect the harmonised data
            variant.accession = f"{variant.accession}_hmPOS_{str(target_build)}"
            yield variant
    else:
        for variant in variants:
            # can't remap, so don't try
            yield variant


def check_effect_allele(variants, drop_missing=False):
    """
    Odd effect allele:

    >>> from pgscatalog.core.lib import models
    >>> variant = models.ScoreVariant(**{"effect_allele": "Z", "effect_weight": 5, "accession": "test", "row_nr": 0, "chr_name": "1", "chr_position": 1})
    >>> list(check_effect_allele([variant], drop_missing=True)) # doctest: +ELLIPSIS
    []
    >>> variant = models.ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0, "chr_name": "1", "chr_position": 1})
    >>> list(check_effect_allele([variant], drop_missing=True)) # doctest: +ELLIPSIS
    [ScoreVariant(..., effect_allele=Allele(allele='A', is_snp=True), ...)]
    """
    n_bad = 0
    for variant in variants:
        if not variant.effect_allele.is_snp:
            n_bad += 1
            if drop_missing:
                continue

        yield variant

    if n_bad > 1:
        logger.warning(f"{n_bad} variants have invalid effect alleles (not ACTG)")


def detect_complex(variants):
    """Some older scoring files in the PGS Catalog are complicated.
    They often require bespoke set up to support interaction terms, etc
    This function only exists to provide loud warnings to end users.
    """
    is_complex = False

    for variant in variants:
        if not is_complex:
            if variant.is_complex:
                is_complex = True

        yield variant

    if is_complex:
        logger.warning("Complex scoring file detected")
        logger.warning(
            "Complex files are difficult to calculate properly and may require manual intervention"
        )


def lift(
    *, scoring_file, harmonised, current_build, target_build, chain_dir, min_lift=0.95
):
    variants = scoring_file.variants

    skip_lo = True
    if target_build != current_build:
        if not harmonised:
            skip_lo = False

    if skip_lo:
        logger.info("Skipping liftover")
        for variant in variants:
            yield variant
    else:
        logger.info("Starting liftover")
        lo = load_chain(
            current_build=current_build, target_build=target_build, chain_dir=chain_dir
        )

        n_lifted = 0
        n = 0

        for variant in variants:
            chrom = "chr" + variant.chr_name  # ew
            pos = int(variant.chr_position) - 1  # VCF -> 1 based, UCSC -> 0 based
            lifted = lo.convert_coordinate(chrom, pos)
            if lifted:
                variant.chr_name = lifted[0][0][3:].split("_")[0]
                variant.chr_position = lifted[0][1] + 1  # reverse 0 indexing
                yield variant
                n_lifted += 1
            else:
                variant.chr_name = None
                variant.chr_position = None
                yield variant
            n += 1

        if (n_lifted / n) < min_lift:
            raise LiftoverError(f"{scoring_file!r}")
        else:
            logger.info("Liftover successful")


def load_chain(*, current_build, target_build, chain_dir):
    """Only supports loading GRCh37 and GRCh38 chain files

    >>> from pgscatalog.core.lib import Config
    >>> chain_dir = Config.ROOT_DIR / "tests" / "data" / "chain"
    >>> load_chain(current_build=GenomeBuild.GRCh37, target_build=GenomeBuild.GRCh38, chain_dir=chain_dir) # doctest: +ELLIPSIS
    <pyliftover.liftover.LiftOver object at...

    >>> load_chain(current_build=GenomeBuild.GRCh38, target_build=GenomeBuild.GRCh37, chain_dir=chain_dir) # doctest: +ELLIPSIS
    <pyliftover.liftover.LiftOver object at...

    >>> load_chain(current_build=GenomeBuild.NCBI36, target_build=GenomeBuild.GRCh38, chain_dir=chain_dir)
    Traceback (most recent call last):
    ...
    ValueError: Unsupported liftover current_build=GenomeBuild.NCBI36, target_build=GenomeBuild.GRCh38
    """
    chain_path = pathlib.Path(chain_dir)

    match (current_build, target_build):
        case GenomeBuild.GRCh37, GenomeBuild.GRCh38:
            chain_path = chain_path / "hg19ToHg38.over.chain.gz"
        case GenomeBuild.GRCh38, GenomeBuild.GRCh37:
            chain_path = chain_path / "hg38ToHg19.over.chain.gz"
        case _:
            raise ValueError(f"Unsupported liftover {current_build=}, {target_build=}")

    if not chain_path.exists():
        raise FileNotFoundError

    return pyliftover.LiftOver(str(chain_path))
