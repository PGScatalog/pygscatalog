"""A declared sex, and what the calls say about it (ADR 0014)

A **Declared sex** is QC input only. It never derives ploidy, never fills a ploidy a
format could not express, and never reaches a dosage or a score - ADR 0007's rejection
of a pipeline-wide sex column as a *scoring* input stands in full. Everything here
compares a declaration against ploidy that was observed per call, and stops when the two
disagree.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from pgscatalog.calc.lib.cache.genotypecache import GenotypeCache
from pgscatalog.calc.lib.config import config
from pgscatalog.calc.lib.constants import (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    NO_ALLELE_SENTINEL_VALUE,
)

if TYPE_CHECKING:
    import pathlib

logger = logging.getLogger(__name__)

# plink and UK Biobank's coding. 0 is a legal value everywhere: it is how a user
# declares an unknown sex, and how they exempt a sample from the check (ADR 0014)
UNKNOWN, MALE, FEMALE = 0, 1, 2
VALID_SEX_CODES = frozenset({UNKNOWN, MALE, FEMALE})

SEX_MAP_COLUMNS = ("sample_id", "sex")


def read_sex_map(path: pathlib.Path) -> dict[str, int]:
    """Read a two-column sex map: sample_id, sex.

    Tab separated, with a header naming both columns. Codes are 0 unknown, 1 male,
    2 female. Any other value is refused naming the file, the line and the value - a
    mandatory declaration is worth nothing if the parser guesses at it (ADR 0014).

    ponytail: hand-rolled over open() rather than duckdb read_csv, which is how this
    repo reads its other tabular input. Two columns do not need a query engine, and
    duckdb's parse errors would have to be remapped to name the offending line anyway.
    """
    codes: dict[str, int] = {}

    with path.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        if tuple(header) != SEX_MAP_COLUMNS:
            raise ValueError(
                f"{path}: sex map header must be "
                f"{list(SEX_MAP_COLUMNS)}, found {header}"
            )

        # line 1 was the header, so data starts at 2 and the number in an error message
        # is the number the user's editor shows
        for line_number, line in enumerate(f, start=2):
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(SEX_MAP_COLUMNS):
                raise ValueError(
                    f"{path} line {line_number}: expected "
                    f"{len(SEX_MAP_COLUMNS)} tab separated columns, found {len(fields)}"
                )

            sample_id, code = fields[0], fields[1]

            try:
                sex = int(code)
            except ValueError:
                raise ValueError(
                    f"{path} line {line_number}: sex must be one of "
                    f"{sorted(VALID_SEX_CODES)}, found {code!r}"
                ) from None

            if sex not in VALID_SEX_CODES:
                raise ValueError(
                    f"{path} line {line_number}: sex must be one of "
                    f"{sorted(VALID_SEX_CODES)}, found {code!r}"
                )

            if sample_id in codes:
                raise ValueError(
                    f"{path} line {line_number}: sample {sample_id!r} is declared twice"
                )

            codes[sample_id] = sex

    if not codes:
        raise ValueError(f"{path}: sex map declares no samples")

    return codes


def check_sex_map_covers(
    sex_map: dict[str, int], samples: list[str], source: str
) -> None:
    """Every genotyped sample must carry a declaration.

    Coverage is asymmetric, and the asymmetry is forced by sample batching (ADR 0003).
    A genotyped sample missing from the map is a hard error, because a mandatory
    declaration that is silently incomplete is not a declaration, and those samples
    would otherwise pass the check unexamined with nothing recorded. A sample in the
    map that this shard does not hold is fine and expected: batches see a subset of the
    cohort, so one cohort-wide map is the normal way to run this, and rejecting extra
    rows would make the map a function of how the cohort was chopped up.
    """
    missing = [sample for sample in samples if sample not in sex_map]

    if missing:
        shown = ", ".join(missing[:10])
        suffix = "" if len(missing) <= 10 else f" (and {len(missing) - 10} more)"
        raise ValueError(
            f"{len(missing)} of {len(samples)} samples in {source} are missing from "
            f"the sex map: {shown}{suffix}. Every genotyped sample needs a "
            f"declaration. A map that misses all of them is usually keyed on a "
            f"different sample ID space."
        )


# hemizygosity is deterministic rather than estimated, so this floor is tolerance to bad
# calls, not statistical power: at one informative variant a single miswritten no-call
# inverts the fraction and stops a run over a sample that is fine. Ten with a majority
# rule tolerates four (ADR 0014)
MIN_INFORMATIVE_VARIANTS = 10

CONSISTENT, MISMATCH, EXEMPT, NOT_CONCLUDED = (
    "consistent",
    "mismatch",
    "exempt",
    "not_concluded",
)


@dataclass(frozen=True)
class SampleResult:
    """What the calls said about one sample's declaration"""

    sample_id: str
    declared: int
    hemizygous_fraction: float | None
    n_informative: int
    outcome: str


class SexCheckMismatch(ValueError):
    """The calls contradict the declaration for at least one sample"""


def _shard_fractions(genotypes, samples: list[str]) -> tuple[list[float | None], int]:
    """Per sample, the fraction of its called genotypes that are hemizygous.

    A variant is **informative** when at least one call in the shard is hemizygous.
    A variant diploid in every sample is excluded: it is either pseudoautosomal or a
    cohort with no hemizygous samples in it, and nothing tells those apart without the
    build-specific PAR coordinates ADR 0007 declined to ship.

    This exclusion is not a refinement. A scoring file whose X variants are mostly
    pseudoautosomal makes a man diploid at most of them, landing his fraction near 0 and
    reading him as female - without the rule the check is wrong, not weak.

    Informativeness counts a hemizygous no-call ([255, 254]): the file declared ploidy
    one there, which is the signal being read. The fraction's denominator does not, so
    that call contributes to neither side of the sample's own ratio.
    """
    n_samples = len(samples)
    hemizygous_called = np.zeros(n_samples, dtype=np.int64)
    n_called = np.zeros(n_samples, dtype=np.int64)
    n_informative = 0

    # read in variant-axis slices rather than whole: a biobank chrX shard is
    # (n_variants, n_samples, 2) uint8, and a PGS can carry thousands of X positions
    for start in range(0, genotypes.shape[0], config.ZARR_VARIANT_CHUNK_SIZE):
        gt = genotypes[start : start + config.ZARR_VARIANT_CHUNK_SIZE]

        # ploidy is derived from the genotype array, never stored (ADR 0007)
        ploidy = (gt != NO_ALLELE_SENTINEL_VALUE).sum(axis=2)
        called = ~(gt == MISSING_GENOTYPE_SENTINEL_VALUE).any(axis=2)
        hemizygous = ploidy == 1

        informative = hemizygous.any(axis=1)
        n_informative += int(informative.sum())
        hemizygous_called += (hemizygous & called)[informative].sum(axis=0)
        n_called += called[informative].sum(axis=0)

    # a NULL fraction, never 0.0: computed naively 0/0 is 0.0, which reads as diploid,
    # which reads as female-consistent - every woman would pass silently and every man
    # be flagged, on a check that never ran (ADR 0014)
    fractions: list[float | None] = [
        float(hemizygous_called[i]) / float(n_called[i]) if n_called[i] else None
        for i in range(n_samples)
    ]

    return fractions, n_informative


def check_sex(
    cache_dir: pathlib.Path, sampleset: str, sex_map: dict[str, int]
) -> list[SampleResult]:
    """Compare every declaration against the ploidy the chrX calls actually carry.

    Runs per shard, with no combine: every (sample, variant) pair appears in at most
    one shard and multiple files under one chromosome are sample batches with disjoint
    samples (ADR 0003), so one chrX shard holds every cached X call for its samples.
    """
    cache = GenotypeCache(cache_dir / "genotypes.zarr", sampleset=sampleset)
    shards = [shard for shard in cache.shards if shard.chromosome == "X"]

    if not shards:
        # the map was still mandatory and goes unread: the mandate is decided from the
        # scoring files, before anything knows which chromosomes will be loaded
        logger.warning(
            "No chrX among the targets, so the sex check has nothing to run on. "
            "The sex map was required because a scoring file carries X or Y variants, "
            "and it has not been used."
        )
        return []

    root = cache.root
    results: list[SampleResult] = []

    for shard in shards:
        samples = shard.samples
        fractions, n_informative = _shard_fractions(
            root[f"{shard.path}/genotypes"], samples
        )

        if n_informative < MIN_INFORMATIVE_VARIANTS:
            logger.warning(
                f"{shard.path}: {n_informative} informative X variants, fewer than "
                f"{MIN_INFORMATIVE_VARIANTS}, so the sex check does not conclude for "
                f"any of its {len(samples)} samples. Not concluding is not a pass."
            )

        for sample_id, fraction in zip(samples, fractions, strict=True):
            declared = sex_map[sample_id]
            results.append(
                SampleResult(
                    sample_id=sample_id,
                    declared=declared,
                    hemizygous_fraction=fraction,
                    n_informative=n_informative,
                    outcome=_outcome(declared, fraction, n_informative),
                )
            )

    return results


def _outcome(declared: int, fraction: float | None, n_informative: int) -> str:
    """One sample's verdict.

    Classified on the majority: declared male expects 1, declared female expects 0.
    A fraction of exactly 0.5 is a mismatch rather than a pass - half the calls say
    male and half say female, so the sample or the file is broken either way, and that
    is contradictory evidence rather than the absent evidence the floor covers.
    """
    if declared == UNKNOWN:
        return EXEMPT

    if n_informative < MIN_INFORMATIVE_VARIANTS or fraction is None:
        return NOT_CONCLUDED

    if fraction == 0.5:
        return MISMATCH

    observed = MALE if fraction > 0.5 else FEMALE
    return CONSISTENT if observed == declared else MISMATCH


def write_sex_check_report(path: pathlib.Path, results: list[SampleResult]) -> None:
    """The per-sample table, written whether or not the run goes on to raise"""
    with path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["sample_id", "declared", "hemizygous_fraction", "n_informative", "outcome"]
        )
        for result in results:
            writer.writerow(
                [
                    result.sample_id,
                    result.declared,
                    ""
                    if result.hemizygous_fraction is None
                    else result.hemizygous_fraction,
                    result.n_informative,
                    result.outcome,
                ]
            )

    logger.info(f"Sex check report written to {path}")


def report_sex_check(results: list[SampleResult]) -> None:
    """Say how much was checked, then raise if anything disagreed.

    The count is unconditional. A user facing four hundred hard failures will set the
    whole column to 0 and re-run, and the pipeline then exits cleanly having verified
    almost nothing - worse than report-only, because the exit code now asserts
    something. A run that checked 3 of 500,000 samples has to say so in the same breath
    as it says it passed (ADR 0014).
    """
    checked = [x for x in results if x.outcome in (CONSISTENT, MISMATCH)]
    exempted = [x for x in results if x.outcome == EXEMPT]
    mismatches = [x for x in results if x.outcome == MISMATCH]

    logger.warning(
        f"Sex check: {len(checked)} of {len(results)} samples checked, "
        f"{len(exempted)} exempted by a declaration of 0."
    )

    if mismatches:
        # named all at once: failing on the first would make a cohort with four hundred
        # mismatches four hundred round trips (ADR 0014)
        named = ", ".join(
            f"{x.sample_id} (declared {x.declared}, hemizygous fraction "
            f"{x.hemizygous_fraction:.2f})"
            for x in mismatches
        )
        raise SexCheckMismatch(
            f"{len(mismatches)} samples whose chrX calls contradict their declaration: "
            f"{named}. Correct the sex map, or declare a sample 0 to exempt it, then "
            f"re-run with --zarr_zip_file to reuse the cache this run already wrote."
        )
