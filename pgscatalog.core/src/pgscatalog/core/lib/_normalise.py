""" This module contains a data processing pipeline to format ScoreVariants in a
standard way. Each step in the data processing pipeline is a generator that operates
on a list of ScoreVariants and yields updated ScoreVariants. This makes it easy to
plug in extra steps where needed, and lazily works on millions of objects."""

import logging
import pathlib

import pyliftover

from .genomebuild import GenomeBuild
from .scorevariant import EffectType, ScoreVariant, EffectAllele
from .pgsexceptions import LiftoverError

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
            harmonised=scoring_file.harmonised,
            current_build=scoring_file.genome_build,
            target_build=target_build,
            chain_dir=chain_dir,
        )
    else:
        variants = scoring_file.variants

    variants = remap_harmonised(variants, scoring_file.harmonised)
    variants = check_bad_variant(variants, drop_missing)

    if drop_missing:
        variants = drop_hla(variants)

    variants = assign_effect_type(variants)
    variants = check_effect_weight(variants)
    variants = assign_other_allele(variants)
    variants = check_effect_allele(variants, drop_missing)
    variants = detect_complex(variants)

    if scoring_file.is_wide:
        # wide data must be sorted because check_duplicates requires sorted input
        variants = (x for x in sorted(variants, key=lambda x: x["accession"]))

    variants = check_duplicates(variants)

    return variants


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

        # None other allele -> empty string
        variant_id: str = ":".join(
            [
                str(getattr(variant, k) or "")
                for k in ["chr_name", "chr_position", "effect_allele", "other_allele"]
            ]
        )

        if variant_id in seen_ids:
            variant.is_duplicated = True
            n_duplicates += 1

        seen_ids[variant_id] = True

        yield variant
        n_variants += 1

    if n_duplicates > 0:
        logger.warning(
            f"{n_duplicates} of {n_variants} variants are duplicated in: {current_accession}"
        )


def drop_hla(variants):
    """Drop HLA alleles from a list of ScoreVariants

    >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0})
    >>> list(drop_hla([variant])) # doctest: +ELLIPSIS
    [ScoreVariant(effect_allele='A',...)]

    >>> variant = ScoreVariant(**{"effect_allele": "P", "effect_weight": 5, "accession": "test", "row_nr": 0})
    >>> list(drop_hla([variant]))
    []
    """
    n_dropped = 0
    for variant in variants:
        match variant:
            case _ if variant.effect_allele in (EffectAllele("P"), EffectAllele("N")):
                n_dropped += 1
                continue
            case _:
                yield variant

    logger.warning(f"{n_dropped} HLA alleles detected and dropped")


def check_effect_weight(variants):
    """Check that effect weights are valid floats. Effect weights are intentionally
    left as strings during processing.

    >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0})
    >>> list(check_effect_weight([variant])) # doctest: +ELLIPSIS
    [ScoreVariant(effect_allele='A',effect_weight=5,...)]

    >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": "potato", "accession": "test", "row_nr": 0})
    >>> list(check_effect_weight([variant])) # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    ValueError
    """
    for variant in variants:
        try:
            float(variant.effect_weight)
        except ValueError as e:
            logger.critical(f"{variant} has bad effect weight")
            raise ValueError from e
        else:
            yield variant


def assign_other_allele(variants):
    """Check if there's more than one possible other allele, remove if true

    >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0, "other_allele": "A"})
    >>> list(assign_other_allele([variant])) # doctest: +ELLIPSIS
    [ScoreVariant(effect_allele='A',...,other_allele='A',...)]
    >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0, "other_allele": "A/C"})
    >>> list(assign_other_allele([variant])) # doctest: +ELLIPSIS
    [ScoreVariant(effect_allele='A',...,other_allele=None,...)]
    """
    n_dropped = 0
    for variant in variants:
        if "/" in variant.other_allele:
            n_dropped += 1
            variant.other_allele = None

        yield variant

    if n_dropped > 0:
        logger.warning(f"Multiple other_alleles detected in {n_dropped} variants")
        logger.warning("Other allele for these variants is set to missing")


def assign_effect_type(variants):
    """Convert PGS Catalog effect type columns to EffectType enums

    The most common type of effect type is additive:

    >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0, "is_recessive": "False", "is_dominant": "False"})
    >>> list(assign_effect_type([variant])) # doctest: +ELLIPSIS
    [ScoreVariant(...,effect_type=EffectType.ADDITIVE,...)]
    >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0, "is_recessive": "True", "is_dominant": "False"})
    >>> list(assign_effect_type([variant])) # doctest: +ELLIPSIS
    [ScoreVariant(...,effect_type=EffectType.RECESSIVE,...)]
    >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0, "is_recessive": "False", "is_dominant": "True"})
    >>> list(assign_effect_type([variant])) # doctest: +ELLIPSIS
    [ScoreVariant(...,effect_type=EffectType.DOMINANT,...)]

    is_recessive and is_dominant fields are parsed from strings to bools during __init__.
    """
    for variant in variants:
        match (variant.is_recessive, variant.is_dominant):
            case (None, None) | (False, False):
                pass  # default value is additive, pass to break match and yield
            case (False, True):
                variant.effect_type = EffectType.DOMINANT
            case (True, False):
                variant.effect_type = EffectType.RECESSIVE
            case _:
                logger.critical(f"Bad effect type setting: {variant}")
                raise Exception
        yield variant


def remap_harmonised(variants, harmonised):
    """
    Overwrite key attributes with harmonised data, if available.

    In this case chr_name, chr_position, and other allele are missing.
    Perhaps authors submitted rsID and effect allele originally:

    >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0, "hm_chr": 1, "hm_pos": 100, "hm_inferOtherAllele": "A"})
    >>> list(remap_harmonised([variant], harmonised=True)) # doctest: +ELLIPSIS
    [ScoreVariant(...,chr_name=1,chr_position=100,...other_allele='A'...)]
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
            yield variant
    else:
        for variant in variants:
            # can't remap, so don't try
            yield variant


def check_bad_variant(variants, drop_missing=False):
    """
    Missing effect allele:

    >>> variant = ScoreVariant(**{"effect_allele": None, "effect_weight": 5, "accession": "test", "row_nr": 0})
    >>> list(check_bad_variant([variant], drop_missing=True)) # doctest: +ELLIPSIS
    []

    Missing chromosome name and position:

    >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0})
    >>> list(check_bad_variant([variant], drop_missing=True)) # doctest: +ELLIPSIS
    []
    """
    n_bad = 0
    for variant in variants:
        match variant:
            case (
                ScoreVariant(chr_name=None)
                | ScoreVariant(chr_position=None)
                | ScoreVariant(effect_allele=None)
            ):
                # (effect weight checked separately)
                n_bad += 1
                if not drop_missing:
                    yield variant
            case _:
                yield variant

    if n_bad > 1:
        logger.warning(f"{n_bad} bad variants")


def check_effect_allele(variants, drop_missing=False):
    """
    Odd effect allele:

    >>> variant = ScoreVariant(**{"effect_allele": "Z", "effect_weight": 5, "accession": "test", "row_nr": 0})
    >>> list(check_effect_allele([variant], drop_missing=True)) # doctest: +ELLIPSIS
    []
    >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0})
    >>> list(check_effect_allele([variant], drop_missing=True)) # doctest: +ELLIPSIS
    [ScoreVariant(effect_allele='A'...)]
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

    >>> from ._config import Config
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
