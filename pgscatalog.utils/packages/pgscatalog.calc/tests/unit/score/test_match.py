import duckdb
import pytest

from pgscatalog.calc.lib.score._matchvariants import (
    ExclusionReason,
    MatchPriority,
    MatchResult,
    _setup_enums,
    check_nucleotide_allele,
    classify_match,
    create_allele_match_table,
    exclude_match_collisions,
    get_variant_match_priority,
    match_variant,
    update_match_table,
)
from pgscatalog.calc.lib.scorefile import load_scoring_files

# Test functions of: update_match_table, _setup_enums, match_variant,
# get_variant_match_priority, check_nucleotide_allele,
# classify_match (complement, complement_table)


@pytest.mark.parametrize(
    "allele",
    ["N", "X", "AN"],
    ids=[
        "N raises ValueError",
        "X raises ValueError",
        "AN raises ValueError",
    ],
)
def test_check_nucleotide_allele(allele: str) -> None:
    """Unmatchable alleles are skipped before this, so what's left is corrupt"""
    with pytest.raises(ValueError):
        check_nucleotide_allele(allele)


@pytest.mark.parametrize("allele", ["A", "GATC"], ids=["one base", "many bases"])
def test_check_nucleotide_allele_accepts_sequences(allele: str) -> None:
    assert check_nucleotide_allele(allele) is None


@pytest.mark.parametrize(
    "effect_allele, other_allele, target_ref, target_alt, result",
    [
        # result: (MatchPriority, effect allele index, matched as written,
        # matched complemented)
        ("A", "C", "A", "C", (MatchPriority.REFALT, 0, True, False)),
        ("A", "C", "C", "A", (MatchPriority.ALTREF, 1, True, False)),
        ("A", "C", "T", "G", (MatchPriority.REFALT_FLIP, 0, False, True)),
        ("A", "C", "G", "T", (MatchPriority.ALTREF_FLIP, 1, False, True)),
        ("A", None, "A", "T", (MatchPriority.REF_NO_OA, 0, True, True)),
        ("A", None, "T", "A", (MatchPriority.ALT_NO_OA, 1, True, True)),
        ("A", None, "C", "T", (MatchPriority.NO_MATCH, None, False, True)),
        ("A", "T", "A", "T", (MatchPriority.REFALT, 0, True, True)),
        ("A", "C", "C", "G", (MatchPriority.NO_MATCH, None, False, False)),
        # docs/adr/0011: the other strand carries the *reverse* complement
        ("AC", "CC", "GT", "GG", (MatchPriority.REFALT_FLIP, 0, False, True)),
        ("AC", "CC", "GG", "GT", (MatchPriority.ALTREF_FLIP, 1, False, True)),
        ("AC", "CC", "TG", "GG", (MatchPriority.NO_MATCH, None, False, False)),
        ("A", "AT", "A", "AT", (MatchPriority.REFALT, 0, True, False)),
        ("A", "AT", "T", "TA", (MatchPriority.NO_MATCH, None, False, False)),
        ("AT", None, "A", "AT", (MatchPriority.ALT_NO_OA, 1, True, False)),
    ],
    ids=[
        "match ref alt",
        "match alt ref",
        "match ref alt flipped",
        "match alt ref flipped",
        "match ref missing other allele",
        "match alt missing other allele",
        # no flipped no-OA match type exists, so this reading only counts towards
        # ambiguity - and here nothing else at the site matched, so it stays a
        # non-match rather than becoming one
        "complemented match missing other allele",
        # both readings match: the score variant as written wins the priority, and
        # the caller is told a second reading matched
        "palindrome matches both ways round",
        "no match",
        # AC read on the other strand is GT, not TG: a plain complement matches the
        # sequences a strand flip doesn't produce and misses the ones it does
        "a multi-base strand flip matches the reference allele",
        "a multi-base strand flip matches the alternative allele",
        "a plain complement of a multi-base target doesn't match",
        # an indel's anchor base is positional, so the strand-flipped record sits
        # somewhere else entirely and no complemented reading is attempted
        "an indel matches forward",
        "an indel never matches on a complemented reading",
        "an indel that is its own reverse complement still isn't read complemented",
    ],
)
def test_classify_match(
    effect_allele: str,
    other_allele: str | None,
    target_ref: str,
    target_alt: str,
    result: tuple,
    alt_allele_idx: int = 1,
) -> None:
    """Test standard matching scenarios"""
    assert (
        classify_match(
            effect_allele=effect_allele,
            other_allele=other_allele,
            target_ref=target_ref,
            target_alt=target_alt,
            alt_allele_index=alt_allele_idx,
        )
        == result
    )


def _best_match(
    effect_allele: str,
    other_allele: str | None,
    ref: str,
    alts: list[str],
    match_ambiguous: bool = False,
    match_flipped: bool = True,
) -> MatchResult:
    """The match a score variant gets at a target variant, as match_variant picks it"""
    return max(
        get_variant_match_priority(
            effect_allele=effect_allele,
            other_allele=other_allele,
            ref=ref,
            alts=alts,
            match_ambiguous=match_ambiguous,
            match_flipped=match_flipped,
        ),
        key=lambda match: match.match_priority.value,
    )


@pytest.mark.parametrize(
    "ref, alts, effect_allele, other_allele, match_type, effect_allele_idx",
    [
        ("A", ["T"], "A", "T", MatchPriority.EXCLUDED, None),
        ("C", ["G"], "C", "G", MatchPriority.EXCLUDED, None),
        ("A", ["C"], "A", "C", MatchPriority.REFALT, 0),
        ("A", ["C"], "C", None, MatchPriority.ALT_NO_OA, 1),
        ("A", ["C", "G"], "C", None, MatchPriority.EXCLUDED, None),
        ("A", ["C", "T"], "A", None, MatchPriority.EXCLUDED, None),
        ("A", ["T", "C"], "C", "A", MatchPriority.ALTREF, 2),
        ("A", ["C"], "T", "G", MatchPriority.REFALT_FLIP, 0),
        ("A", ["T"], "G", "C", MatchPriority.NO_MATCH, None),
        ("A", ["TA"], "A", "TA", MatchPriority.REFALT, 0),
        ("AT", ["TA"], "AT", "TA", MatchPriority.EXCLUDED, None),
        ("GT", ["GG"], "AC", "CC", MatchPriority.REFALT_FLIP, 0),
        ("TG", ["GG"], "AC", "CC", MatchPriority.NO_MATCH, None),
        ("T", ["TA"], "A", "AT", MatchPriority.NO_MATCH, None),
        ("A", ["AT"], "AT", None, MatchPriority.ALT_NO_OA, 1),
    ],
    ids=[
        # docs/adr/0010: ambiguity is a property of the pairing of a score variant
        # with a target variant, never of a pair of alleles
        "A/T at a palindromic target is ambiguous, as before",
        "C/G at a palindromic target is ambiguous, as before",
        "A/C is resolvable",
        "a no-OA score variant is resolvable when only one reading matches",
        "a no-OA score variant is ambiguous when the target carries the complement",
        "a no-OA reference match is ambiguous when the target carries the complement",
        "a palindromic pair at the target leaves other score variants resolvable",
        "a resolved strand flip is not ambiguity: only one reading matched",
        # this one used to be EXCLUDED, because ambiguity was computed from the
        # target's alleles before anyone asked whether the score variant matched
        "a score variant that matches nothing at a palindromic target is unmatched",
        "A/TA is resolvable",
        # AT and TA are each their own reverse complement, so this pair really is
        # palindromic - it stays ambiguous under docs/adr/0011, for the right reason
        # rather than the accidental one
        "AT/TA matches both ways round, so it's ambiguous",
        # docs/adr/0011: revcomp(AC) = GT, not TG
        "a multi-base strand flip resolves",
        "a plain complement of a multi-base target isn't a strand flip",
        "an indel is never read complemented, so it can't flip-match",
        # under a naive reverse complement this would be ambiguous: revcomp(AT) = AT
        # is the alternative allele the score variant already matched forward
        "an indel is never read complemented, so it's never ambiguous",
    ],
)
def test_ambiguity_is_a_property_of_the_pairing(
    ref: str,
    alts: list[str],
    effect_allele: str,
    other_allele: str | None,
    match_type: MatchPriority,
    effect_allele_idx: int | None,
) -> None:
    match = _best_match(effect_allele, other_allele, ref, alts)

    assert match.match_priority == match_type
    assert match.effect_allele_idx == effect_allele_idx
    assert match.is_ambiguous == (match_type == MatchPriority.EXCLUDED)


@pytest.mark.parametrize(
    "effect_allele, other_allele, match_type, effect_allele_idx",
    [
        ("A", "C", MatchPriority.REFALT, 0),
        ("A", "G", MatchPriority.REFALT, 0),
        ("C", "A", MatchPriority.ALTREF, 1),
        ("G", "A", MatchPriority.ALTREF, 2),
    ],
    ids=[
        "effect A other C matches the reference allele",
        "effect A other G matches the reference allele",
        "effect C other A matches alternative allele 1",
        "effect G other A matches alternative allele 2",
    ],
)
def test_every_alternative_allele_is_a_match_candidate(
    effect_allele: str,
    other_allele: str,
    match_type: MatchPriority,
    effect_allele_idx: int,
) -> None:
    """REF=A ALTS=[C,G]: matching used to stop after the first alternative allele,
    so the two score variants pairing with G were reported unmatched."""
    match = _best_match(effect_allele, other_allele, "A", ["C", "G"])

    assert match.match_priority == match_type
    assert match.effect_allele_idx == effect_allele_idx
    assert match.is_multiallelic is True


def test_alternative_allele_order_doesnt_change_the_answer() -> None:
    """ALT column order carries no meaning, so neither may the result.

    A/T at REF=A ALTS=[C,T] used to report unmatched, because alt C yielded
    NO_MATCH and matching then stopped. Written ALTS=[T,C] the same target reported
    excluded, because the ambiguous branch skipped past the break.
    """
    assert _best_match("A", "T", "A", ["C", "T"]) == _best_match(
        "A", "T", "A", ["T", "C"]
    )
    assert _best_match("A", "T", "A", ["C", "T"]).match_priority is (
        MatchPriority.EXCLUDED
    )


def test_unmatchable_allele_keeps_its_index() -> None:
    """Alleles are positional, so * is skipped where it stands (docs/adr/0010)"""
    match = _best_match("G", "A", "A", ["*", "G"])

    assert match.match_priority is MatchPriority.ALTREF
    assert match.effect_allele_idx == 2


def test_an_indel_is_never_excluded_as_ambiguous() -> None:
    """The whole decision an indel gets, through the UDF the database calls.

    A VCF indel's shared anchor base is positional - it stays leftmost on whichever
    strand the record is read - so the strand-flipped record sits at a different
    position with different alleles, and no string operation on these two recovers
    it. No complemented reading is attempted, so the score variant matches forward
    and there is no second reading to be ambiguous against (docs/adr/0011).

    Reverse-complementing the indel anyway would exclude this: revcomp(AT) is AT, the
    alternative allele it already matched as written.
    """
    result = match_variant(
        False,  # match_ambiguous
        False,  # match_multiallelic
        True,  # match_flipped
        "AT",  # effect_allele
        None,  # other_allele
        "A",  # ref
        ["AT"],
    )

    assert result["match_type"] == "ALT_NO_OA"
    assert result["effect_allele_idx"] == 1
    assert result["is_ambiguous"] is False
    assert result["exclusion_reason"] is None


@pytest.mark.parametrize(
    "effect_allele, other_allele, ref, alts, match_ambiguous, match_results",
    [
        (
            "A",
            "T",
            "A",
            ["T"],
            True,
            [
                MatchResult(
                    MatchPriority.REFALT,
                    effect_allele_idx=0,
                    is_ambiguous=True,
                    is_multiallelic=False,
                )
            ],
        ),
        (
            "T",
            "A",
            "A",
            ["T"],
            True,
            [
                MatchResult(
                    MatchPriority.ALTREF,
                    effect_allele_idx=1,
                    is_ambiguous=True,
                    is_multiallelic=False,
                )
            ],
        ),
        (
            "A",
            "C",
            "A",
            ["C", "T"],
            False,
            [
                MatchResult(
                    MatchPriority.REFALT,
                    effect_allele_idx=0,
                    is_ambiguous=False,
                    is_multiallelic=True,
                ),
                MatchResult(
                    MatchPriority.NO_MATCH,
                    effect_allele_idx=None,
                    is_ambiguous=False,
                    is_multiallelic=True,
                ),
            ],
        ),
        (
            "A",
            "T",
            "A",
            ["C", "T"],
            False,
            [
                MatchResult(
                    MatchPriority.EXCLUDED,
                    effect_allele_idx=None,
                    is_ambiguous=True,
                    is_multiallelic=True,
                    exclusion_reason=ExclusionReason.AMBIGUOUS,
                )
            ],
        ),
        (
            "A",
            "C",
            "A",
            ["C"],
            False,
            [
                MatchResult(
                    MatchPriority.REFALT,
                    effect_allele_idx=0,
                    is_ambiguous=False,
                    is_multiallelic=False,
                )
            ],
        ),
        (
            "C",
            "A",
            "A",
            ["C"],
            False,
            [
                MatchResult(
                    MatchPriority.ALTREF,
                    effect_allele_idx=1,
                    is_ambiguous=False,
                    is_multiallelic=False,
                )
            ],
        ),
    ],
    ids=[
        "match ambiguous refalt",
        "match ambiguous altref",
        "match multiallelic refalt",
        "excluded ambiguous multiallelic",
        "match refalt",
        "match altref",
    ],
)
def test_match_priority(
    effect_allele: str,
    other_allele: str | None,
    ref: str,
    alts: list[str],
    match_ambiguous: bool,
    match_results: list[MatchResult],
) -> None:
    """Test more complicated matching scenarios.

    One result per alternative allele, so a multiallelic variant yields several and
    the caller keeps the best. An excluded score variant is the exception: the
    decision is about the pairing rather than about one allele, so it's made once.
    """
    result: list[MatchResult] = list(
        get_variant_match_priority(
            effect_allele=effect_allele,
            other_allele=other_allele,
            ref=ref,
            alts=alts,
            match_ambiguous=match_ambiguous,
            match_flipped=True,
        )
    )

    assert result == match_results


@pytest.mark.parametrize(
    "match_dict, is_matched, match_type, match_priority",
    [
        (
            {
                "effect_allele_idx": 0,
                "match_priority": MatchPriority.REFALT,
                "is_ambiguous": False,
                "is_multiallelic": False,
            },
            True,
            "REFALT",
            7,
        ),
        (
            {
                "effect_allele_idx": 1,
                "match_priority": MatchPriority.ALTREF,
                "is_ambiguous": False,
                "is_multiallelic": False,
            },
            True,
            "ALTREF",
            6,
        ),
        (
            {
                "effect_allele_idx": 0,
                "match_priority": MatchPriority.REFALT_FLIP,
                "is_ambiguous": True,
                "is_multiallelic": False,
            },
            True,
            "REFALT_FLIP",
            5,
        ),
        (
            {
                "effect_allele_idx": None,
                "match_priority": MatchPriority.EXCLUDED,
                "is_ambiguous": True,
                "is_multiallelic": False,
                "exclusion_reason": ExclusionReason.AMBIGUOUS,
            },
            False,
            "EXCLUDED",
            1,
        ),
        (
            {
                "effect_allele_idx": None,
                "match_priority": MatchPriority.NO_MATCH,
                "is_ambiguous": False,
                "is_multiallelic": False,
            },
            False,
            "NO_MATCH",
            0,
        ),
    ],
    ids=[
        "refalt matchresult",
        "altref matchresult",
        "refalt_flip matchresult",
        "excluded matchresult",
        "no match matchresult",
    ],
)
def test_matchresult(
    match_dict: dict, is_matched: bool, match_type: str, match_priority: int
):
    """Test computed match result properties (is_matched, match_typpe)"""
    result = MatchResult(**match_dict)
    assert result.is_matched == is_matched
    assert result.match_type == match_type
    assert frozenset(result.to_dict().keys()) == {
        "effect_allele_idx",
        "match_priority",
        "match_type",
        "match_summary",
        "is_ambiguous",
        "is_multiallelic",
        "is_matched",
        "exclusion_reason",
    }
    assert result.to_dict()["match_priority"] == match_priority


@pytest.mark.parametrize(
    "match_dict, result",
    [
        (
            {
                "match_ambiguous": False,
                "match_multiallelic": False,
                "match_flipped": True,
                "effect_allele": "A",
                "other_allele": "G",
                "ref": "A",
                "alts": ["G"],
            },
            {
                "effect_allele_idx": 0,
                "match_priority": 7,
                "match_type": "REFALT",
                "match_summary": "matched",
                "is_ambiguous": False,
                "is_matched": True,
                "is_multiallelic": False,
                "exclusion_reason": None,
            },
        ),
        (
            {
                "match_ambiguous": False,
                "match_multiallelic": False,
                "match_flipped": True,
                "effect_allele": "G",
                "other_allele": "A",
                "ref": "A",
                "alts": ["G"],
            },
            {
                "effect_allele_idx": 1,
                "match_priority": 6,
                "match_type": "ALTREF",
                "match_summary": "matched",
                "is_ambiguous": False,
                "is_matched": True,
                "is_multiallelic": False,
                "exclusion_reason": None,
            },
        ),
        (
            {
                "match_ambiguous": True,
                "match_multiallelic": False,
                "match_flipped": True,
                "effect_allele": "A",
                "other_allele": "T",
                "ref": "A",
                "alts": ["T"],
            },
            {
                "effect_allele_idx": 0,
                "match_priority": 7,
                "match_type": "REFALT",
                "match_summary": "matched",
                "is_ambiguous": True,
                "is_matched": True,
                "is_multiallelic": False,
                "exclusion_reason": None,
            },
        ),
        (
            {
                "match_ambiguous": False,
                "match_multiallelic": False,
                "match_flipped": True,
                "effect_allele": "A",
                "other_allele": "T",
                "ref": "A",
                "alts": ["T"],
            },
            {
                "effect_allele_idx": None,
                "match_priority": 1,
                "match_type": "EXCLUDED",
                "match_summary": "excluded",
                "is_ambiguous": True,
                "is_matched": False,
                "is_multiallelic": False,
                "exclusion_reason": "ambiguous",
            },
        ),
        (
            {
                "match_ambiguous": False,
                "match_multiallelic": False,
                "match_flipped": True,
                "effect_allele": "A",
                "other_allele": "G",
                "ref": "A",
                "alts": ["G", "C"],
            },
            {
                "effect_allele_idx": None,
                "match_priority": 1,
                "match_type": "EXCLUDED",
                "match_summary": "excluded",
                "is_ambiguous": None,
                "is_matched": False,
                "is_multiallelic": True,
                "exclusion_reason": "multiallelic",
            },
        ),
        (
            {
                "match_ambiguous": False,
                "match_multiallelic": True,
                "match_flipped": True,
                "effect_allele": "A",
                "other_allele": "G",
                "ref": "A",
                "alts": ["G", "C"],
            },
            {
                "effect_allele_idx": 0,
                "match_priority": 7,
                "match_type": "REFALT",
                "match_summary": "matched",
                "is_ambiguous": False,
                "is_matched": True,
                "is_multiallelic": True,
                "exclusion_reason": None,
            },
        ),
        (
            {
                "match_ambiguous": False,
                "match_multiallelic": False,
                "match_flipped": True,
                "effect_allele": "A",
                "other_allele": "T",
                "ref": "C",
                "alts": ["C"],
            },
            {
                "effect_allele_idx": None,
                "match_priority": 0,
                "match_type": "NO_MATCH",
                "match_summary": "unmatched",
                "is_ambiguous": False,
                "is_matched": False,
                "is_multiallelic": False,
                "exclusion_reason": None,
            },
        ),
    ],
    ids=[
        "simple_REFALT_match",
        "simple_ALTREF_match",
        "ambiguous_but_allowed",
        "ambiguous_but_disallowed",
        "multiallelic_but_disallowed",
        "multiallelic_but_allowed",
        "no_match",
    ],
)
def test_match_variant_udf(match_dict: dict, result: dict) -> None:
    """Test match_variant, a user defined function for duckdb"""
    assert match_variant(**match_dict) == result


@pytest.mark.parametrize(
    "alts",
    [["*"], ["<NON_REF>"]],
    ids=[
        "spanning deletion is unmatched, not an error",
        "symbolic allele is unmatched, not an error",
    ],
)
def test_match_variant_unmatchable_allele(alts: list[str]) -> None:
    """Legal VCF alleles that name no sequence yield NO_MATCH instead of raising"""
    assert match_variant(
        False,  # match_ambiguous
        False,  # match_multiallelic
        True,  # match_flipped
        "A",  # effect_allele
        "T",  # other_allele
        "A",  # ref
        alts,
    ) == {
        "effect_allele_idx": None,
        "match_priority": 0,
        "match_type": "NO_MATCH",
        "match_summary": "unmatched",
        "is_ambiguous": False,
        "is_matched": False,
        "is_multiallelic": False,
        "exclusion_reason": None,
    }


def test_match_variant_corrupt_allele() -> None:
    """A corrupt allele string isn't unmatchable, it's wrong, and still raises"""
    with pytest.raises(ValueError):
        match_variant(False, False, True, "A", "T", "A", ["XYZ"])


def test_ignoring_strand_flips_excludes_a_flipped_match() -> None:
    """The variant was found, so it's excluded and says why, not reported unmatched.

    revcomp(T/G) is A/C, so this score variant matches REF=A ALT=C only when read
    complemented. Nothing matches it as written (docs/adr/0011).
    """
    kept = _best_match(effect_allele="T", other_allele="G", ref="A", alts=["C"])
    assert kept.match_type == "REFALT_FLIP"

    dropped = _best_match(
        effect_allele="T", other_allele="G", ref="A", alts=["C"], match_flipped=False
    )
    assert dropped.match_summary == "excluded"
    assert dropped.exclusion_reason is ExclusionReason.STRAND_FLIP
    assert dropped.effect_allele_idx is None


def test_ignoring_strand_flips_leaves_a_forward_match_alone() -> None:
    """Only a variant that matched *by* a complemented reading is dropped"""
    assert (
        _best_match(
            effect_allele="A",
            other_allele="C",
            ref="A",
            alts=["C"],
            match_flipped=False,
        ).match_type
        == "REFALT"
    )


def test_ignoring_strand_flips_keeps_the_ambiguity_check() -> None:
    """The reading is still made, so a palindrome is still ambiguous.

    Skipping the complemented reading instead would make is_ambiguous unreachable and
    quietly score every palindromic SNP on a forward-strand assumption, which is the
    opposite of what the flag promises (docs/adr/0011).
    """
    result = _best_match(
        effect_allele="A", other_allele="T", ref="A", alts=["T"], match_flipped=False
    )
    assert result.is_ambiguous is True
    assert result.exclusion_reason is ExclusionReason.AMBIGUOUS


def test_a_kept_ambiguous_palindrome_still_matches_forward() -> None:
    """Both flags raised: the pairing matched both ways round, and the forward
    reading is the one that won, so --ignore_strand_flips has nothing to drop."""
    result = _best_match(
        effect_allele="A",
        other_allele="T",
        ref="A",
        alts=["T"],
        match_ambiguous=True,
        match_flipped=False,
    )
    assert result.is_ambiguous is True
    assert result.match_type == "REFALT"


def test_multiallelic_takes_precedence_over_ambiguous() -> None:
    """Both reasons apply and only one column is available to say so.

    REF=A ALTS=[T,C] carries a palindromic pair, so the score variant A/T is
    ambiguous there as well as multiallelic. match_variant tests multiallelic before
    comparing any allele, so that's the reason reported (docs/adr/0010).
    """
    result = match_variant(
        False,  # match_ambiguous
        False,  # match_multiallelic
        True,  # match_flipped
        "A",  # effect_allele
        "T",  # other_allele
        "A",  # ref
        ["T", "C"],
    )

    assert result["match_summary"] == "excluded"
    assert result["exclusion_reason"] == "multiallelic"


def test_matchresult_excluded_without_a_reason_raises() -> None:
    """An excluded row a reader can't interpret is the bug this column fixes"""
    with pytest.raises(ValueError, match="exclusion_reason"):
        MatchResult(
            MatchPriority.EXCLUDED,
            effect_allele_idx=None,
            is_ambiguous=True,
            is_multiallelic=False,
        )


def test_matchresult_reason_without_an_exclusion_raises() -> None:
    """A reason on a matched row would say a decision was taken that wasn't"""
    with pytest.raises(ValueError, match="exclusion_reason"):
        MatchResult(
            MatchPriority.REFALT,
            effect_allele_idx=0,
            is_ambiguous=False,
            is_multiallelic=False,
            exclusion_reason=ExclusionReason.AMBIGUOUS,
        )


def test_setup_enums() -> None:
    conn = duckdb.connect(":memory:")
    _setup_enums(conn)

    conn = duckdb.connect(":memory:")
    _setup_enums(conn)

    # Check that match_summary_enum exists and has expected values
    (summary_values,) = conn.execute(
        "SELECT enum_range(NULL::match_summary_enum)"
    ).fetchone()
    assert summary_values == ["matched", "unmatched", "excluded"]

    # Check that match_type_enum has the same labels as MatchPriority
    (type_values,) = conn.execute("SELECT enum_range(NULL::match_type_enum)").fetchone()
    assert type_values == [e.name for e in MatchPriority]

    # every exclusion names its reason, and the reasons are enumerable
    (reason_values,) = conn.execute(
        "SELECT enum_range(NULL::exclusion_reason_enum)"
    ).fetchone()
    assert reason_values == [
        "multiallelic",
        "ambiguous",
        "strand_flip",
        "collision",
    ]


def _make_test_conn(
    pgs000001, alts: list[str] | None = None
) -> duckdb.DuckDBPyConnection:
    # nothing here needs the database on disk
    conn = duckdb.connect(":memory:")

    # score_variant_table loading via load_scoring_files functions from Scorefiles
    load_scoring_files(conn=conn, scorefile_paths=[pgs000001])

    # minimal targetvariants table matching your SQL in update_match_table
    conn.execute("""
        CREATE TABLE targetvariants (
            geno_index UINTEGER,
            chr_name TEXT,
            chr_pos UINTEGER,
            ref TEXT,
            alts TEXT[],
            zarr_group TEXT,
            sampleset TEXT
        );
    """)

    # One matching target row in the same sampleset
    conn.execute(
        """
        INSERT INTO targetvariants
        (geno_index, chr_name, chr_pos, ref, alts, zarr_group, sampleset)
        VALUES (?, ?, ?, ?, ?, ?, ?);
        """,
        [
            10,
            "11",
            69516650,
            "T",
            alts if alts is not None else ["C"],
            "test/ALL/file.bgen",
            "test",
        ],
    )

    return conn


def test_update_match_table_spanning_deletion(pgs000001) -> None:
    """A target record holding a spanning deletion at a scoring file position.

    The raise used to happen inside the match_variant UDF, so it surfaced as an
    opaque duckdb query error. The record is now reported unmatched instead.
    """
    conn = _make_test_conn(pgs000001, alts=["*"])

    update_match_table(
        conn=conn,
        match_ambiguous=False,
        match_multiallelic=False,
        match_flipped=True,
        sampleset="test",
    )

    rows = conn.execute(
        """
        SELECT
            match_result.effect_allele_idx AS effect_allele_idx,
            match_result.match_type        AS match_type,
            match_result.match_summary     AS match_summary,
            match_result.is_matched        AS is_matched
        FROM allele_match_table;
        """
    ).fetchall()

    assert rows == [(None, "NO_MATCH", "unmatched", False)]


def test_update_match_table(pgs000001) -> None:
    conn = _make_test_conn(pgs000001)

    # This calls _setup_enums and registers match_variant UDF internally
    update_match_table(
        conn=conn,
        match_ambiguous=False,
        match_multiallelic=False,
        match_flipped=True,
        sampleset="test",
    )

    # Unpack the struct match_result columns to simple scalars
    rows = conn.execute(
        """
        SELECT
            sampleset,
            accession,
            row_nr,
            match_result.effect_allele_idx AS effect_allele_idx,
            match_result.match_priority    AS match_priority,
            match_result.match_type        AS match_type,
            match_result.match_summary     AS match_summary,
            match_result.is_ambiguous      AS is_ambiguous,
            match_result.is_matched        AS is_matched,
            match_result.is_multiallelic   AS is_multiallelic,
            target_row_nr,
            zarr_group
        FROM allele_match_table;
        """
    ).fetchall()

    # We expect exactly one row with a REFALT match
    assert rows == [
        (
            "test",  # sampleset
            "PGS000001",  # accession
            0,  # row_nr
            0,  # effect_allele_idx
            7,  # match_priority (REFALT)
            "REFALT",  # match_type
            "matched",  # match_summary
            False,  # is_ambiguous
            True,  # is_matched
            False,  # is_multiallelic
            10,  # target_row_nr (geno_index)
            "test/ALL/file.bgen",  # zarr_group
        )
    ]


# --- match collisions (docs/adr/0010)
#
# Two score variants of one accession resolving to the same allele of the same target
# variant. Both carry a weight for one dosage, so both are excluded: picking one is
# guessing what the score's author meant.

COLLISION_SHARD = "test/X/x.vcf.gz"

# a scoring file processed by pgscatalog-format, as load_scoring_files loads it
_SCORE_ROW_COLUMNS = 9


def _make_collision_conn(
    score_rows: list[tuple], target_rows: list[tuple]
) -> duckdb.DuckDBPyConnection:
    """A connection holding hand-written score and target variants.

    Collisions need two scoring file rows at one position, which no test scoring
    file carries, so these are written directly into the tables update_match_table
    reads.
    """
    conn = duckdb.connect(":memory:")
    conn.execute("""
        CREATE TABLE score_variant_table (
            chr_name VARCHAR,
            chr_position UINTEGER,
            effect_allele VARCHAR,
            other_allele VARCHAR,
            effect_weight DOUBLE,
            effect_type VARCHAR,
            is_duplicated BOOLEAN,
            accession VARCHAR,
            row_nr UINTEGER,
            PRIMARY KEY(row_nr, accession)
        );
    """)
    placeholders = ", ".join("?" * _SCORE_ROW_COLUMNS)
    conn.executemany(
        f"INSERT INTO score_variant_table VALUES ({placeholders});", score_rows
    )
    conn.execute("""
        CREATE TABLE targetvariants (
            geno_index UINTEGER,
            chr_name TEXT,
            chr_pos UINTEGER,
            ref TEXT,
            alts TEXT[],
            zarr_group TEXT,
            sampleset TEXT
        );
    """)
    conn.executemany(
        "INSERT INTO targetvariants VALUES (?, ?, ?, ?, ?, ?, ?);", target_rows
    )
    return conn


# a scoring file with a genuine XY and X variant at one position. Chromosome names
# are canonicalised on the way into score_variant_table (ADR 0008), so XY has already
# become X here: two rows, one canonical position, two different weights, and
# is_duplicated False on both because the flag was computed upstream on the allele
# pair and these are not duplicates by that definition
XY_AND_X_SCORE_ROWS = [
    ("X", 100, "A", "C", 0.1, "additive", False, "PGS000001", 0),
    ("X", 100, "A", "C", 0.2, "additive", False, "PGS000001", 1),
]
XY_AND_X_TARGET_ROWS = [(5, "X", 100, "A", ["C"], COLLISION_SHARD, "test")]


def _xy_and_x_conn(score_rows: list[tuple] | None = None) -> duckdb.DuckDBPyConnection:
    return _make_collision_conn(
        score_rows=XY_AND_X_SCORE_ROWS if score_rows is None else score_rows,
        target_rows=XY_AND_X_TARGET_ROWS,
    )


def _match_outcomes(conn: duckdb.DuckDBPyConnection) -> list[tuple]:
    return conn.execute("""
        SELECT accession,
               row_nr,
               match_result.effect_allele_idx AS effect_allele_idx,
               match_result.match_type        AS match_type,
               match_result.match_summary     AS match_summary,
               match_result.is_matched        AS is_matched
        FROM allele_match_table
        ORDER BY accession, row_nr;
    """).fetchall()


def test_collision_excludes_both_score_variants() -> None:
    """Neither weight is chosen, because choosing is guessing.

    Both rows used to survive matching, and the weight matrix pivot then kept one of
    them with any_value - an unspecified choice between two weights for one dosage.
    """
    conn = _xy_and_x_conn()

    update_match_table(
        conn=conn,
        match_ambiguous=False,
        match_multiallelic=False,
        match_flipped=True,
        sampleset="test",
    )

    assert _match_outcomes(conn) == [
        ("PGS000001", 0, None, "EXCLUDED", "excluded", False),
        ("PGS000001", 1, None, "EXCLUDED", "excluded", False),
    ]

    # with no effect allele index there's no weight to pivot: create_wide_weights_table
    # drops both rows, so neither weight can reach the score
    (weighted,) = conn.execute("""
        SELECT count(*) FROM allele_match_table
        WHERE match_result.effect_allele_idx IS NOT NULL;
    """).fetchone()
    assert weighted == 0


def test_collision_outcome_is_the_same_on_every_run() -> None:
    """The bug was that it wasn't: both rows tied on every ordering the pipeline
    applied, so which weight survived was up to the query engine."""
    outcomes = []
    # the two rows are loaded in the opposite order the second time round: nothing
    # about the outcome may depend on which one the engine happens to see first
    for score_rows in (XY_AND_X_SCORE_ROWS, list(reversed(XY_AND_X_SCORE_ROWS))):
        conn = _xy_and_x_conn(score_rows)
        update_match_table(
            conn=conn,
            match_ambiguous=False,
            match_multiallelic=False,
            match_flipped=True,
            sampleset="test",
        )
        outcomes.append(_match_outcomes(conn))

    assert outcomes[0] == outcomes[1]
    assert {row[4] for row in outcomes[0]} == {"excluded"}


def test_two_score_variants_at_one_position_are_not_always_a_collision() -> None:
    """One score variant matches the reference allele, the other the alternative.

    Two alleles of one variant, two variants found: they take different dosages, so
    there's nothing to guess and both are scored.
    """
    conn = _make_collision_conn(
        score_rows=[
            ("X", 100, "A", "C", 0.1, "additive", False, "PGS000001", 0),
            ("X", 100, "C", "A", 0.2, "additive", False, "PGS000001", 1),
        ],
        target_rows=[(5, "X", 100, "A", ["C"], COLLISION_SHARD, "test")],
    )

    update_match_table(
        conn=conn,
        match_ambiguous=False,
        match_multiallelic=False,
        match_flipped=True,
        sampleset="test",
    )

    assert _match_outcomes(conn) == [
        ("PGS000001", 0, 0, "REFALT", "matched", True),
        ("PGS000001", 1, 1, "ALTREF", "matched", True),
    ]


def test_match_rate_counts_a_target_allele_once() -> None:
    """The min_overlap numerator counts rows with match_summary 'matched'.

    Nothing dedups the match table - its primary key is per scoring variant - so a
    colliding pair used to put two matched rows on one target allele and count one
    locus twice. Excluding collisions is what removes that: no second mechanism
    needed (ADR 0010). Two score variants on two *different* alleles of one variant
    stay two matched rows, because they are two variants found.
    """
    conn = _make_collision_conn(
        score_rows=[
            # collide on the reference allele of target row 5
            ("X", 100, "A", "C", 0.1, "additive", False, "PGS000001", 0),
            ("X", 100, "A", "C", 0.2, "additive", False, "PGS000001", 1),
            # a clean match at another locus
            ("X", 200, "G", "T", 0.3, "additive", False, "PGS000001", 2),
        ],
        target_rows=[
            (5, "X", 100, "A", ["C"], COLLISION_SHARD, "test"),
            (6, "X", 200, "G", ["T"], COLLISION_SHARD, "test"),
        ],
    )

    update_match_table(
        conn=conn,
        match_ambiguous=False,
        match_multiallelic=False,
        match_flipped=True,
        sampleset="test",
    )

    numerator = conn.execute("""
        SELECT target_row_nr, count(*)
        FROM allele_match_table
        WHERE match_result.match_summary = 'matched'
        GROUP BY accession, zarr_group, target_row_nr,
                 match_result.effect_allele_idx
        ORDER BY target_row_nr;
    """).fetchall()

    # the colliding locus contributes nothing, and the locus that did match counts
    # exactly once
    assert numerator == [(6, 1)]


def _match_table_with_rows(rows: list[tuple]) -> duckdb.DuckDBPyConnection:
    """A match table filled directly, to exercise the collision pass on its own.

    (accession, row_nr, effect_allele_idx, target_row_nr, zarr_group)
    """
    conn = duckdb.connect(":memory:")
    _setup_enums(conn)
    create_allele_match_table(conn, target_alleles=True)
    conn.executemany(
        """
        INSERT INTO allele_match_table VALUES (
            'test',
            ?,
            ?,
            struct_pack(
                effect_allele_idx := ?::UINT8,
                is_multiallelic := FALSE,
                is_ambiguous := FALSE,
                is_matched := TRUE,
                match_priority := 7::UINT8,
                match_type := 'REFALT'::match_type_enum,
                match_summary := 'matched'::match_summary_enum,
                exclusion_reason := NULL::exclusion_reason_enum
            ),
            ?,
            ?,
            'A',
            ['C']
        );
        """,
        rows,
    )
    return conn


def test_different_alternative_alleles_are_not_a_collision() -> None:
    """Two score variants matching alternative alleles 1 and 2 of one multiallelic
    variant. Two alleles of one variant, so two variants found: they take different
    dosages and nothing has to be guessed (ADR 0010)."""
    conn = _match_table_with_rows(
        [
            ("PGS000001", 0, 1, 7, COLLISION_SHARD),
            ("PGS000001", 1, 2, 7, COLLISION_SHARD),
        ]
    )

    assert exclude_match_collisions(conn) == 0
    assert [row[2] for row in _match_outcomes(conn)] == [1, 2]


def test_one_target_row_nr_in_two_shards_is_two_variants() -> None:
    """target_row_nr is a shard-local geno_index, so the same number in two shards
    names two different variants and can't collide."""
    conn = _match_table_with_rows(
        [
            ("PGS000001", 0, 0, 7, COLLISION_SHARD),
            ("PGS000001", 1, 0, 7, "test/X/other.vcf.gz"),
        ]
    )

    assert exclude_match_collisions(conn) == 0
    assert [row[5] for row in _match_outcomes(conn)] == [True, True]


def test_two_accessions_on_one_target_allele_are_not_a_collision() -> None:
    """Each accession gets its own weight matrix column, so there's no competition"""
    conn = _match_table_with_rows(
        [
            ("PGS000001", 0, 0, 7, COLLISION_SHARD),
            ("PGS000002", 0, 0, 7, COLLISION_SHARD),
        ]
    )

    assert exclude_match_collisions(conn) == 0
    assert [row[5] for row in _match_outcomes(conn)] == [True, True]


def test_collision_keeps_the_target_variant_flags() -> None:
    """is_ambiguous and is_multiallelic describe the target variant, not the
    decision, so exclusion leaves them alone and names itself in exclusion_reason."""
    conn = _match_table_with_rows(
        [
            ("PGS000001", 0, 0, 7, COLLISION_SHARD),
            ("PGS000001", 1, 0, 7, COLLISION_SHARD),
        ]
    )
    conn.execute("""
        UPDATE allele_match_table
        SET match_result = struct_pack(
            effect_allele_idx := match_result.effect_allele_idx,
            is_multiallelic := TRUE,
            is_ambiguous := FALSE,
            is_matched := match_result.is_matched,
            match_priority := match_result.match_priority,
            match_type := match_result.match_type,
            match_summary := match_result.match_summary,
            exclusion_reason := match_result.exclusion_reason
        );
    """)

    assert exclude_match_collisions(conn) == 2

    flags = conn.execute("""
        SELECT match_result.is_multiallelic, match_result.is_ambiguous,
               match_result.match_priority, match_result.exclusion_reason
        FROM allele_match_table;
    """).fetchall()
    # a multiallelic target that collided reports the collision, not the flag: the
    # flag was already true when the row matched
    assert flags == [(True, False, MatchPriority.EXCLUDED.value, "collision")] * 2


# every reason, in one match table, written by both writers: match_variant sets the
# first two from inside the UDF, the collision pass sets the third across rows

_EVERY_REASON_SCORE_ROWS = [
    # two score variants resolving to the reference allele of target row 5
    ("X", 100, "A", "C", 0.1, "additive", False, "PGS000001", 0),
    ("X", 100, "A", "C", 0.2, "additive", False, "PGS000001", 1),
    # a palindromic pairing: both readings match target row 6
    ("X", 200, "A", "T", 0.3, "additive", False, "PGS000001", 2),
    # target row 7 carries two alternative alleles
    ("X", 300, "A", "C", 0.4, "additive", False, "PGS000001", 3),
    # a clean match, and a score variant that matches nothing at target row 9
    ("X", 400, "G", "T", 0.5, "additive", False, "PGS000001", 4),
    ("X", 500, "A", "G", 0.6, "additive", False, "PGS000001", 5),
]
_EVERY_REASON_TARGET_ROWS = [
    (5, "X", 100, "A", ["C"], COLLISION_SHARD, "test"),
    (6, "X", 200, "A", ["T"], COLLISION_SHARD, "test"),
    (7, "X", 300, "A", ["C", "G"], COLLISION_SHARD, "test"),
    (8, "X", 400, "G", ["T"], COLLISION_SHARD, "test"),
    (9, "X", 500, "G", ["T"], COLLISION_SHARD, "test"),
]


def _reasons(conn: duckdb.DuckDBPyConnection) -> list[tuple]:
    return conn.execute("""
        SELECT row_nr,
               match_result.match_summary    AS match_summary,
               match_result.exclusion_reason AS exclusion_reason
        FROM allele_match_table
        ORDER BY row_nr;
    """).fetchall()


def test_every_exclusion_names_its_reason() -> None:
    """All three reasons at once: excluded rows name one, everything else is NULL"""
    conn = _make_collision_conn(
        score_rows=_EVERY_REASON_SCORE_ROWS,
        target_rows=_EVERY_REASON_TARGET_ROWS,
    )

    update_match_table(
        conn=conn,
        match_ambiguous=False,
        match_multiallelic=False,
        match_flipped=True,
        sampleset="test",
    )

    assert _reasons(conn) == [
        (0, "excluded", "collision"),
        (1, "excluded", "collision"),
        (2, "excluded", "ambiguous"),
        (3, "excluded", "multiallelic"),
        (4, "matched", None),
        (5, "unmatched", None),
    ]


def test_keeping_multiallelics_leaves_the_other_reasons_alone() -> None:
    """--keep_multiallelic empties one reason and doesn't touch the other two"""
    conn = _make_collision_conn(
        score_rows=_EVERY_REASON_SCORE_ROWS,
        target_rows=_EVERY_REASON_TARGET_ROWS,
    )

    update_match_table(
        conn=conn,
        match_ambiguous=False,
        match_multiallelic=True,
        match_flipped=True,
        sampleset="test",
    )

    assert _reasons(conn) == [
        (0, "excluded", "collision"),
        (1, "excluded", "collision"),
        (2, "excluded", "ambiguous"),
        # the multiallelic target is matched now, and still flagged multiallelic
        (3, "matched", None),
        (4, "matched", None),
        (5, "unmatched", None),
    ]
    assert conn.execute("""
        SELECT match_result.is_multiallelic FROM allele_match_table WHERE row_nr = 3;
    """).fetchone() == (True,)


def test_unmatched_rows_never_collide() -> None:
    """Every unmatched row at one locus has a NULL effect allele index. They aren't
    two weights competing for one dosage - they're no weights at all."""
    conn = _match_table_with_rows(
        [
            ("PGS000001", 0, None, 7, COLLISION_SHARD),
            ("PGS000001", 1, None, 7, COLLISION_SHARD),
        ]
    )

    assert exclude_match_collisions(conn) == 0
