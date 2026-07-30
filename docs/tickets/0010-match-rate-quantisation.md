# 0010 — Match rate quantises before summing, and can flip the gate

**Refs:** 0004

## Why

`_matchlog.py` casts each subgroup's `fraction` to `DECIMAL(5, 2)` *before* summing
them into the match rate:

```sql
CAST(count / sum(count) OVER (PARTITION BY sampleset, accession) AS DECIMAL(5, 2))
    AS fraction
```

The `matched` rows are split across several subgroups (by `is_ambiguous` and
`is_multiallelic`), so each subgroup's fraction is rounded to two decimal places and
the roundings accumulate. Three matched subgroups at 245/1000 each sum to `0.75` when
the true rate is `0.735`. The error is up to ±0.02 and can move the rate either way.

Combined with the strict `>` comparison against `min_overlap`, an accession sitting
exactly on the threshold fails.

## Scope

Sum the raw fractions and round once, at the end — or better, divide the summed counts.
Rounding is a presentation concern for the summary log, not an input to the gate.

## Out of scope

Not fixed alongside the alpha2 split, because it changes match log output and would
have invalidated 0004's "match logs match today's output" acceptance criterion.

## Acceptance

- An accession whose true match rate is 0.735 across three subgroups is rejected at
  `min_overlap = 0.75`.
- The `match_rate` column still displays at two decimal places.
