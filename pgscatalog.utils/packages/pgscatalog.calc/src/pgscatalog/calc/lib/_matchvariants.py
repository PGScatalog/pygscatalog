from __future__ import annotations

import enum
import functools
import itertools
import logging
from dataclasses import dataclass
from functools import cache
from typing import TYPE_CHECKING, ClassVar

import duckdb
from duckdb.sqltypes import BOOLEAN, UINTEGER, VARCHAR, DuckDBPyType

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from collections.abc import Callable, Generator
    from typing import Literal


def _setup_enums(conn: duckdb.DuckDBPyConnection) -> None:
    conn.execute("""
        CREATE TYPE IF NOT EXISTS match_summary_enum AS
        ENUM ('matched', 'unmatched', 'excluded');
    """)
    conn.execute(f"""
        CREATE TYPE IF NOT EXISTS match_type_enum AS
        ENUM{tuple(e.name for e in MatchPriority)}
    """)


def update_match_table(
    conn: duckdb.DuckDBPyConnection,
    match_ambiguous: bool,
    match_multiallelic: bool,
    sampleset: str,
) -> None:
    logger.info("Matching score variants with target variants by position and alleles")
    # python UDF function parameters have to get values from db columns
    # use functools.partial to "pre-fill" the function with non-column parameters
    # be careful to use posargs not kwargs (duckdb complains about kwargs)
    match_variant_function: Callable[
        [bool, bool, str, str | None, str, list[str]], dict
    ] = functools.partial(match_variant, match_ambiguous, match_multiallelic)
    # register a Python user defined function
    conn.create_function(
        "match_variant",
        match_variant_function,
        return_type=MatchResult.struct_type,
    )

    _setup_enums(conn=conn)
    logger.info("Creating temporary table with position matches")
    conn.execute(
        """
        CREATE TEMP TABLE position_match_table AS
        SELECT
            score.accession,
            score.row_nr,
            score.effect_allele,
            score.other_allele,
            target.geno_index AS target_row_nr,
            target.chr_name,
            target.chr_pos,
            target.ref,
            target.alts,
            target.filename,
            target.sampleset
        FROM score_variant_table AS score
        JOIN target_variants.variants_table AS target
        ON score.chr_name = target.chr_name
        AND score.chr_position = target.chr_pos
        WHERE target.ref IS NOT NULL AND target.sampleset = ?;
        """,
        [sampleset],
    )

    logger.info("Creating allele match table if it doesn't exist")
    conn.sql("""
        -- custom allele matching logic is applied here
        CREATE TABLE IF NOT EXISTS allele_match_table (
            sampleset TEXT NOT NULL,
            accession TEXT NOT NULL,
            row_nr UINTEGER NOT NULL,
            match_result STRUCT(effect_allele_idx UINT8,
                is_multiallelic BOOLEAN,
                is_ambiguous BOOLEAN,
                is_matched BOOLEAN,
                match_priority UINT8,
                match_type match_type_enum,
                match_summary match_summary_enum) NOT NULL,
            target_row_nr UINTEGER NOT NULL,
            filename TEXT NOT NULL,
            PRIMARY KEY (sampleset, accession, row_nr)
        );
        """)

    logger.info("Matching alleles to find effect allele index")
    conn.sql("""
        -- first apply the allele matching UDF
        WITH match_candidates AS (
            SELECT
                sampleset,
                accession,
                row_nr,
                filename,
                match_variant(effect_allele, other_allele, ref, alts)
                    AS match_result,
                target_row_nr
            FROM position_match_table
        ),
        -- filter rows during insert to only include the best match_priority
        -- breaks ties with row_number
        ranked_matches AS (
            SELECT *,
                   ROW_NUMBER() OVER (
                       PARTITION BY accession, row_nr
                       ORDER BY match_result.match_priority DESC
                   ) AS rn
            FROM match_candidates
        )
        INSERT INTO allele_match_table
        SELECT
            sampleset,
            accession,
            row_nr,
            match_result,
            target_row_nr,
            filename
        FROM ranked_matches
        WHERE rn = 1
        ORDER BY accession, row_nr;
        """)
    pass


def match_variant(
    match_ambiguous: bool,
    match_multiallelic: bool,
    effect_allele: str,
    other_allele: str | None,
    ref: str,
    alts: list[str],
) -> dict:
    # effect allele index
    # if effect allele = REF and gt is (0, 0) for a sample then dosage is 2
    best_match: MatchResult
    is_multiallelic: bool = False
    if alts is not None and len(alts) > 1:
        is_multiallelic = True

    if is_multiallelic and not match_multiallelic:
        best_match = MatchResult(
            match_priority=MatchPriority.EXCLUDED,
            effect_allele_idx=None,
            is_ambiguous=None,
            is_multiallelic=is_multiallelic,
        )
    else:
        # in biallelic cases the generator yields one result: the best possible match
        # in multiallelic cases the generator yields one result for each alt allele
        matches = get_variant_match_priority(
            ref=ref,
            alts=alts,
            effect_allele=effect_allele,
            other_allele=other_allele,
            match_ambiguous=match_ambiguous,
        )
        # always only interested in the best possible match for each variant
        best_match = max(matches, key=lambda x: x.match_priority.value)

    return best_match.to_dict()


class MatchPriority(enum.Enum):
    REFALT = 7
    ALTREF = 6
    REFALT_FLIP = 5
    ALTREF_FLIP = 4
    REF_NO_OA = 3
    ALT_NO_OA = 2
    EXCLUDED = 1
    NO_MATCH = 0


@dataclass(slots=True, frozen=True)
class MatchResult:
    struct_type: ClassVar[DuckDBPyType] = duckdb.struct_type(
        {
            "effect_allele_idx": UINTEGER,
            "is_multiallelic": BOOLEAN,
            "is_ambiguous": BOOLEAN,
            "is_matched": BOOLEAN,
            "match_priority": UINTEGER,
            "match_type": VARCHAR,
            "match_summary": VARCHAR,
        }
    )

    match_priority: MatchPriority
    effect_allele_idx: int | None
    is_ambiguous: bool | None
    is_multiallelic: bool | None

    @property
    def is_matched(self) -> bool:
        return self.effect_allele_idx is not None and self.match_priority.value > 1

    @property
    def match_summary(self) -> Literal["matched", "unmatched", "excluded"]:
        match self.match_priority:
            case MatchPriority.EXCLUDED:
                return "excluded"
            case MatchPriority.NO_MATCH:
                return "unmatched"
            case _:
                return "matched"

    @property
    def match_type(self) -> str:
        return self.match_priority.name

    def to_dict(self) -> dict:
        return {
            "effect_allele_idx": self.effect_allele_idx,
            "match_priority": self.match_priority.value,
            "match_type": self.match_type,
            "is_ambiguous": self.is_ambiguous,
            "is_matched": self.is_matched,
            "is_multiallelic": self.is_multiallelic,
            "match_summary": self.match_summary,
        }


def classify_match(
    *,
    effect_allele: str,
    other_allele: str | None,
    target_ref: str,
    target_alt: str | None,
    alt_allele_index: int,
) -> tuple[MatchPriority, int | None]:
    complemented_effect_allele: str = complement(effect_allele)
    complemented_other_allele: str | None = (
        complement(other_allele) if other_allele is not None else None
    )

    if other_allele is not None:
        if effect_allele == target_ref and other_allele == target_alt:
            return MatchPriority.REFALT, 0

        if effect_allele == target_alt and other_allele == target_ref:
            return MatchPriority.ALTREF, alt_allele_index

        if (
            complemented_effect_allele == target_ref
            and complemented_other_allele == target_alt
        ):
            return MatchPriority.REFALT_FLIP, 0

        if (
            complemented_effect_allele == target_alt
            and complemented_other_allele == target_ref
        ):
            return MatchPriority.ALTREF_FLIP, alt_allele_index
    else:
        if effect_allele == target_ref:
            return MatchPriority.REF_NO_OA, 0

        if effect_allele == target_alt:
            return MatchPriority.ALT_NO_OA, alt_allele_index

    return MatchPriority.NO_MATCH, None


def get_variant_match_priority(
    *,
    effect_allele: str,
    other_allele: str | None,
    ref: str,
    alts: list[str],
    match_ambiguous: bool,
) -> Generator[MatchResult, None, None]:
    is_multiallelic: bool = len(alts) > 1

    for alt_allele_idx, (target_ref, target_alt) in enumerate(
        itertools.product([ref], alts), start=1
    ):
        is_ambiguous = is_biallelic_variant_ambiguous(ref=target_ref, alt=target_alt)
        effect_allele_idx: int | None
        priority: MatchPriority

        if is_ambiguous and not match_ambiguous:
            yield MatchResult(
                match_priority=MatchPriority.EXCLUDED,
                effect_allele_idx=None,
                is_ambiguous=is_ambiguous,
                is_multiallelic=is_multiallelic,
            )
            continue

        priority, effect_allele_idx = classify_match(
            effect_allele=effect_allele,
            other_allele=other_allele,
            target_ref=target_ref,
            target_alt=target_alt,
            alt_allele_index=alt_allele_idx,
        )

        yield MatchResult(
            effect_allele_idx=effect_allele_idx,
            is_ambiguous=is_ambiguous,
            match_priority=priority,
            is_multiallelic=is_multiallelic,
        )

        if is_multiallelic:
            # explicitly abandon matching additional alleles for now
            break


@cache
def is_biallelic_variant_ambiguous(ref: str, alt: str) -> bool:
    """
    Identify ambiguous variants (A/T & C/G SNPs)

    Often scores and genotypes will be on different builds and including these SNPs
    would involve tallying incorrect dosages if there has been a strand-flip across
    builds. Even on the same build you can get improperly strand-oriented data.
    """
    valid = frozenset({"A", "C", "T", "G"})

    if not frozenset(ref).issubset(valid):
        raise ValueError(ref)

    if not frozenset(alt).issubset(valid):
        raise ValueError(alt)

    if len(ref) == 1 and len(alt) == 1:
        return complement(ref) == alt

    return False


@cache
def complement_table() -> dict[int, int]:
    return str.maketrans("ATCG", "TAGC")


@cache
def complement(allele: str) -> str:
    return allele.translate(complement_table())
