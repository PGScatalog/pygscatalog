# 0011 — `match_flipped` has always counted zero

**Refs:** 0004

## Why

`_matchlog.py` counts strand-flipped matches with:

```sql
COUNT(*) FILTER (WHERE match_type LIKE '%flip') AS match_flipped
```

`match_type` is `match_type_enum`, whose values are uppercase (`REFALT_FLIP`,
`ALTREF_FLIP`), and duckdb's `LIKE` is case-sensitive. The predicate has never matched
anything, so `match_flipped` reads 0 in every summary log ever produced — including for
data where strand flips genuinely happened. It is exactly the diagnostic a user reaches
for when a score looks wrong because of a build or strand mismatch.

## Scope

`ILIKE` instead of `LIKE`, or compare against the enum values directly.

## Out of scope

Not fixed alongside the alpha2 split, because it changes match log output and would
have invalidated 0004's "match logs match today's output" acceptance criterion.

## Acceptance

- A sampleset with strand-flipped matches reports a non-zero `match_flipped`.
