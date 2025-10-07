import pytest

from pgscatalog.calc.lib._matchvariants import (
    MatchPriority,
    MatchResult,
    classify_match,
    get_variant_match_priority,
    is_biallelic_variant_ambiguous,
    match_variant,
)


@pytest.mark.parametrize(
    "ref, alt, result",
    [
        ("A", "T", True),
        ("C", "G", True),
        ("A", "C", False),
        ("G", "T", False),
        ("AT", "TA", False),
        ("A", "TA", False),
    ],
    ids=[
        "A/T is ambiguous",
        "C/G is ambiguous",
        "A/C is not ambiguous",
        "G/T is not ambiguous",
        "AT/TA is not ambiguous",
        "A/TA is not ambiguous",
    ],
)
def test_match_variant_ambiguous(ref: str, alt: str, result: bool) -> None:
    assert is_biallelic_variant_ambiguous(ref, alt) == result


@pytest.mark.parametrize(
    "ref, alt",
    [("N", "T"), ("X", "A"), ("A", "N")],
    ids=[
        "N/T raises ValueError",
        "X/A raises ValueError",
        "A/N raises ValueError",
    ],
)
def test_bad_ambiguous(ref: str, alt: str) -> None:
    with pytest.raises(ValueError):
        is_biallelic_variant_ambiguous(ref, alt)


@pytest.mark.parametrize(
    "effect_allele, other_allele, target_ref, target_alt, result",
    [
        # result: (MatchPriority, effect allele index which is 0 or 1 for bialleles)
        ("A", "C", "A", "C", (MatchPriority.REFALT, 0)),
        ("A", "C", "C", "A", (MatchPriority.ALTREF, 1)),
        ("A", "C", "T", "G", (MatchPriority.REFALT_FLIP, 0)),
        ("A", "C", "G", "T", (MatchPriority.ALTREF_FLIP, 1)),
        ("A", None, "A", "T", (MatchPriority.REF_NO_OA, 0)),
        ("A", None, "T", "A", (MatchPriority.ALT_NO_OA, 1)),
        ("A", "C", "C", "G", (MatchPriority.NO_MATCH, None)),
    ],
    ids=[
        "match ref alt",
        "match alt ref",
        "match ref alt flipped",
        "match alt ref flipped",
        "match ref missing other allele",
        "match alt missing other allele",
        "no match",
    ],
)
def test_classify_match(
    effect_allele: str,
    other_allele: str | None,
    target_ref: str,
    target_alt: str,
    result: MatchResult,
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


@pytest.mark.parametrize(
    "effect_allele, other_allele, ref, alts, match_ambiguous, match_result",
    [
        (
            "A",
            "T",
            "A",
            ["T"],
            True,
            MatchResult(
                MatchPriority.REFALT,
                effect_allele_idx=0,
                is_ambiguous=True,
                is_multiallelic=False,
            ),
        ),
        (
            "T",
            "A",
            "A",
            ["T"],
            True,
            MatchResult(
                MatchPriority.ALTREF,
                effect_allele_idx=1,
                is_ambiguous=True,
                is_multiallelic=False,
            ),
        ),
        (
            "A",
            "C",
            "A",
            ["C", "T"],
            False,
            MatchResult(
                MatchPriority.REFALT,
                effect_allele_idx=0,
                is_ambiguous=False,
                is_multiallelic=True,
            ),
        ),
        (
            "A",
            "C",
            "A",
            ["C"],
            False,
            MatchResult(
                MatchPriority.REFALT,
                effect_allele_idx=0,
                is_ambiguous=False,
                is_multiallelic=False,
            ),
        ),
        (
            "C",
            "A",
            "A",
            ["C"],
            False,
            MatchResult(
                MatchPriority.ALTREF,
                effect_allele_idx=1,
                is_ambiguous=False,
                is_multiallelic=False,
            ),
        ),
    ],
    ids=[
        "match ambiguous refalt",
        "match ambiguous altref",
        "match multiallelic refalt",
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
    match_result: MatchResult,
) -> None:
    """Test more complicated matching scenarios"""
    result: list[MatchResult] = list(
        get_variant_match_priority(
            effect_allele=effect_allele,
            other_allele=other_allele,
            ref=ref,
            alts=alts,
            match_ambiguous=match_ambiguous,
        )
    )

    # when implemented, multiallelic variants will return one match result for each
    # alternate allele. currently only the first alternate allele is matched
    assert len(result) == 1
    assert result[0] == match_result


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
    }
    assert result.to_dict()["match_priority"] == match_priority


@pytest.mark.parametrize(
    "match_dict, result",
    [
        (
            {
                "match_ambiguous": False,
                "match_multiallelic": False,
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
            },
        ),
        (
            {
                "match_ambiguous": False,
                "match_multiallelic": False,
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
            },
        ),
        (
            {
                "match_ambiguous": True,
                "match_multiallelic": False,
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
            },
        ),
        (
            {
                "match_ambiguous": False,
                "match_multiallelic": False,
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
            },
        ),
        (
            {
                "match_ambiguous": False,
                "match_multiallelic": False,
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
            },
        ),
        (
            {
                "match_ambiguous": False,
                "match_multiallelic": True,
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
            },
        ),
        (
            {
                "match_ambiguous": False,
                "match_multiallelic": False,
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
