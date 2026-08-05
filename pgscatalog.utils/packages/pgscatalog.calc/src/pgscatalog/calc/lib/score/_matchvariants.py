from __future__ import annotations

import enum
import functools
import logging
from dataclasses import dataclass
from functools import cache
from typing import TYPE_CHECKING, ClassVar, NamedTuple

import duckdb
from duckdb.sqltypes import BOOLEAN, UINTEGER, VARCHAR, DuckDBPyType

from pgscatalog.calc.lib.constants import is_unmatchable_allele

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
        # important or missing other alleles will return NULL automatically
        null_handling="special",  # type: ignore
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
            target.zarr_group,
            target.sampleset
        FROM score_variant_table AS score
        JOIN targetvariants AS target
        ON score.chr_name = target.chr_name
        AND score.chr_position = target.chr_pos
        WHERE target.ref IS NOT NULL AND target.sampleset = ?;
        """,
        [sampleset],
    )

    logger.info("Creating allele match table if it doesn't exist")
    create_allele_match_table(conn, target_alleles=True)

    logger.info("Matching alleles to find effect allele index")
    conn.sql("""
        -- first apply the allele matching UDF
        WITH match_candidates AS (
            SELECT
                sampleset,
                accession,
                row_nr,
                zarr_group,
                match_variant(effect_allele, other_allele, ref, alts)
                    AS match_result,
                target_row_nr,
                ref AS target_ref,
                alts AS target_alts
            FROM position_match_table
        ),
        /*
        keep the best match_priority *per leaf*, not per scoring variant. Under
        sample batching the same (accession, row_nr) legitimately matches in every
        batch, and discarding all but one would leave the losing batch's shard with
        no weight-matrix row: its samples would silently vanish from the scores.
        Dedup to one row per scoring variant happens once, in the match log.
        target_row_nr breaks priority ties deterministically.
        */
        ranked_matches AS (
            SELECT *,
                   ROW_NUMBER() OVER (
                       PARTITION BY accession, row_nr, zarr_group
                       ORDER BY match_result.match_priority DESC, target_row_nr
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
            zarr_group,
            target_ref,
            target_alts
        FROM ranked_matches
        WHERE rn = 1
        ORDER BY accession, row_nr, zarr_group;
        """)

    exclude_match_collisions(conn)


def exclude_match_collisions(conn: duckdb.DuckDBPyConnection) -> int:
    """Exclude score variants that resolved to the same target allele.

    A match collision is two score variants of one accession whose best matches land
    on the same allele of the same target variant, so two weights compete for one
    dosage. Choosing between them means guessing what the score's author meant, so
    both are excluded instead (docs/adr/0010). Splitting multiallelics upstream of a
    GWAS, chromosome canonicalisation (docs/adr/0008) and build remapping all produce
    them, and the scoring file looks perfectly well formed: the rows differ in their
    alleles, so `is_duplicated` rightly says they aren't duplicates.

    This runs after matching because a collision isn't a property of a score variant.
    It's only visible across rows, once every score variant has its best match, so
    the earliest place it can be seen is here - which also means the match table this
    stage exports is collision-free everywhere downstream.

    Best-match selection can't do this job, and the two look alike. Best match runs
    many candidate matches -> one score variant, and both mechanisms above (the
    `max()` in the UDF, the ranking window) collapse candidates that way. A collision
    runs the other way, many score variants -> one target allele, and each colliding
    row already *is* the best match for its own score variant.

    The key is (accession, target_row_nr, effect_allele_idx), qualified with the
    sampleset and the shard because target_row_nr is a shard-local geno_index: the
    same number in two shards is two different variants. Two score variants matching
    two *different* alternative alleles of one multiallelic variant differ in
    effect_allele_idx, so they aren't a collision - they're two variants found.

    Returns the number of excluded rows.
    """
    logger.info("Excluding score variants that collide on one target allele")
    conn.execute("""
        CREATE OR REPLACE TEMP TABLE match_collision_table AS
        SELECT sampleset,
               accession,
               zarr_group,
               target_row_nr,
               match_result.effect_allele_idx AS effect_allele_idx
        FROM allele_match_table
        -- an unmatched or already excluded row has no allele to collide on
        WHERE match_result.effect_allele_idx IS NOT NULL
        GROUP BY ALL
        HAVING count(*) > 1;
    """)

    # the struct is rebuilt in the column's declared field order. is_ambiguous and
    # is_multiallelic describe the target variant and stay as matching found them;
    # only the decision changes. When exclusions grow a reason column, this is where
    # a collision names itself.
    excluded: int = conn.execute(
        f"""
        UPDATE allele_match_table AS m
        SET match_result = struct_pack(
            effect_allele_idx := NULL::UINT8,
            is_multiallelic := m.match_result.is_multiallelic,
            is_ambiguous := m.match_result.is_ambiguous,
            is_matched := FALSE,
            match_priority := {MatchPriority.EXCLUDED.value}::UINT8,
            match_type := '{MatchPriority.EXCLUDED.name}'::match_type_enum,
            match_summary := 'excluded'::match_summary_enum
        )
        FROM match_collision_table AS c
        WHERE m.sampleset = c.sampleset
          AND m.accession = c.accession
          AND m.zarr_group = c.zarr_group
          AND m.target_row_nr = c.target_row_nr
          AND m.match_result.effect_allele_idx = c.effect_allele_idx;
        """
    ).fetchone()[0]  # type: ignore[index]

    if excluded:
        logger.warning(
            f"{excluded} score variants excluded: two or more scoring file variants "
            "of one accession matched the same target allele, so their weights "
            "can't be told apart"
        )

    return excluded


def create_allele_match_table(
    conn: duckdb.DuckDBPyConnection, *, target_alleles: bool
) -> None:
    """The match table, declared as the stage that's about to fill it consumes it.

    target_ref / target_alts are carried here rather than joined back from
    targetvariants, which only MATCH ever populates. MATCH writes them and
    MATCH_COMBINE reads them into the match log, so both pass
    ``target_alleles=True``.

    SCORE reads neither, and passes ``target_alleles=False`` on purpose: the
    ungated ``matches.parquet`` carries both columns, so ``INSERT ... BY NAME``
    refuses it against this narrower table.
    """
    target_allele_columns = (
        """
            target_ref TEXT,
            target_alts TEXT[],
        """
        if target_alleles
        else ""
    )
    conn.sql(f"""
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
            zarr_group TEXT NOT NULL,
            {target_allele_columns}
            -- one match per scoring variant per shard (docs/adr/0003)
            PRIMARY KEY (sampleset, accession, row_nr, zarr_group)
        );
        """)


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
        # the generator yields one result for each alternative allele the target
        # carries, so a multiallelic variant yields several. Only one of them can be
        # a real match (docs/adr/0010), so max() picks between a match and the
        # non-matches around it, never between two matches
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


class AlleleMatch(NamedTuple):
    """One score variant read against one alternative allele of a target variant.

    Both readings are reported, not just the winning one, because ambiguity is
    computed from whether a second reading matched (docs/adr/0010) and the priority
    alone can't say: when a score variant matches as written *and* complemented, the
    forward reading is the one that wins.
    """

    priority: MatchPriority
    effect_allele_idx: int | None
    # the score variant as written matches this pair of alleles
    matched_forward: bool
    # a complemented reading of the score variant matches this pair of alleles
    matched_complemented: bool


def _read_alleles(
    *,
    effect_allele: str,
    other_allele: str | None,
    target_ref: str,
    target_alt: str | None,
    alt_allele_index: int,
) -> tuple[MatchPriority, int | None]:
    """Compare one reading of a score variant with one pair of target alleles.

    Sequence equality only: complementing is the caller's job, so that the forward
    and complemented readings stay distinguishable.
    """
    if other_allele is not None:
        if effect_allele == target_ref and other_allele == target_alt:
            return MatchPriority.REFALT, 0

        if effect_allele == target_alt and other_allele == target_ref:
            return MatchPriority.ALTREF, alt_allele_index
    else:
        if effect_allele == target_ref:
            return MatchPriority.REF_NO_OA, 0

        if effect_allele == target_alt:
            return MatchPriority.ALT_NO_OA, alt_allele_index

    return MatchPriority.NO_MATCH, None


# a complemented reading that matches is a resolved strand flip, and only when the
# score variant has an other allele to resolve it with. Without one there's no
# flipped match type to report, so a complemented no-OA reading counts towards
# ambiguity and nothing else
_FLIPPED: dict[MatchPriority, MatchPriority] = {
    MatchPriority.REFALT: MatchPriority.REFALT_FLIP,
    MatchPriority.ALTREF: MatchPriority.ALTREF_FLIP,
}


def classify_match(
    *,
    effect_allele: str,
    other_allele: str | None,
    target_ref: str,
    target_alt: str | None,
    alt_allele_index: int,
) -> AlleleMatch:
    """Read a score variant against one pair of target alleles, both ways round.

    The complemented reading is attempted whether or not the score variant has an
    other allele, which it wasn't before: a no-OA score variant of C matches
    alternative allele 1 of REF=A ALTS=[C,G] as written and alternative allele 2
    complemented, and then nothing says which strand it was written on. That case
    needs the whole allele list to see, so the caller decides what a second matching
    reading means (docs/adr/0010) - this only reports that there was one.
    """
    forward_priority, forward_idx = _read_alleles(
        effect_allele=effect_allele,
        other_allele=other_allele,
        target_ref=target_ref,
        target_alt=target_alt,
        alt_allele_index=alt_allele_index,
    )
    complemented_priority, complemented_idx = _read_alleles(
        effect_allele=complement(effect_allele),
        other_allele=complement(other_allele) if other_allele is not None else None,
        target_ref=target_ref,
        target_alt=target_alt,
        alt_allele_index=alt_allele_index,
    )

    matched_forward: bool = forward_priority is not MatchPriority.NO_MATCH
    matched_complemented: bool = complemented_priority is not MatchPriority.NO_MATCH

    if matched_forward:
        # the score variant as written wins over a complemented reading of it
        return AlleleMatch(
            forward_priority, forward_idx, matched_forward, matched_complemented
        )

    if (flipped := _FLIPPED.get(complemented_priority)) is not None:
        return AlleleMatch(flipped, complemented_idx, matched_forward, True)

    return AlleleMatch(
        MatchPriority.NO_MATCH, None, matched_forward, matched_complemented
    )


def get_variant_match_priority(
    *,
    effect_allele: str,
    other_allele: str | None,
    ref: str,
    alts: list[str],
    match_ambiguous: bool,
) -> Generator[MatchResult, None, None]:
    """Match a score variant against every alternative allele of a target variant.

    One MatchResult is yielded per alternative allele and the caller keeps the best
    priority. No tie-break policy is needed, because at most one alternative allele
    can match a given score variant: with an other allele, a forward match on alt_i
    and a complemented match on alt_j both require alt_i = complement(ref) = alt_j;
    without one, the effect allele equals exactly one distinct alternative allele.
    max() therefore never arbitrates between two real matches (docs/adr/0010).

    Ambiguity is a property of the pairing rather than of a pair of alleles: the
    score variant is ambiguous when a complemented reading of it also matches an
    allele the target carries. That's why every alternative allele is read both ways
    round before anything is yielded - a no-OA score variant can match one
    alternative allele forward and another one complemented, which no per-allele
    check could see. A palindromic pair of alleles at the target is no longer
    ambiguous by itself: it only makes the score variants that match it ambiguous.
    """
    is_multiallelic: bool = len(alts) > 1
    matches: list[AlleleMatch] = []

    for alt_allele_idx, target_alt in enumerate(alts, start=1):
        if is_unmatchable_allele(target_alt):
            # legal VCF, but there's no sequence here to compare a score variant
            # with, so no comparison and no ambiguity check: the allele keeps its
            # index in the variant and yields NO_MATCH (docs/adr/0010)
            matches.append(
                AlleleMatch(
                    priority=MatchPriority.NO_MATCH,
                    effect_allele_idx=None,
                    matched_forward=False,
                    matched_complemented=False,
                )
            )
            continue

        # unmatchable alleles are skipped above, so anything left that isn't a
        # nucleotide sequence is a corrupt allele string, and comparing it would
        # quietly report a non-match
        check_nucleotide_allele(ref)
        check_nucleotide_allele(target_alt)

        matches.append(
            classify_match(
                effect_allele=effect_allele,
                other_allele=other_allele,
                target_ref=ref,
                target_alt=target_alt,
                alt_allele_index=alt_allele_idx,
            )
        )

    # both readings land on an allele the target carries, so the score variant
    # doesn't say which strand it was written on
    is_ambiguous: bool = any(match.matched_forward for match in matches) and any(
        match.matched_complemented for match in matches
    )

    if is_ambiguous and not match_ambiguous:
        # the decision is about the score variant, not about one of the alleles it
        # was compared with, so it's made once
        yield MatchResult(
            match_priority=MatchPriority.EXCLUDED,
            effect_allele_idx=None,
            is_ambiguous=True,
            is_multiallelic=is_multiallelic,
        )
        return

    for match in matches:
        yield MatchResult(
            effect_allele_idx=match.effect_allele_idx,
            is_ambiguous=is_ambiguous,
            match_priority=match.priority,
            is_multiallelic=is_multiallelic,
        )


@cache
def check_nucleotide_allele(allele: str) -> None:
    """Raise ValueError on an allele that isn't a nucleotide sequence.

    Callers must skip unmatchable alleles first (see is_unmatchable_allele): those
    are legal VCF and are kept in place, so what's left to raise on is a corrupt
    allele string.
    """
    valid = frozenset({"A", "C", "T", "G"})

    if not frozenset(allele).issubset(valid):
        raise ValueError(allele)


@cache
def complement_table() -> dict[int, int]:
    return str.maketrans("ATCG", "TAGC")


@cache
def complement(allele: str) -> str:
    return allele.translate(complement_table())
