# 0031 — `EXCLUDED` means three different things and the log doesn't say which

**Refs:** ADR 0010.

## Why

`match_summary = 'excluded'` covers three unrelated decisions, and a reader has to infer
which one happened from `is_ambiguous` and `is_multiallelic`. Those flags describe the
*variant*, not the *decision*, and they are not interchangeable: `is_multiallelic = true`
on a matched row is a fact about the target, on an excluded row it is the cause. After
ticket 0027 there is a third reason with no flag at all, and the inference stops working
entirely.

The summary log is where this bites. `temp_summary_log_table` groups on
`match_summary, is_ambiguous, is_multiallelic` (`lib/score/_matchlog.py:114-118`), so a
user is told "217 variants excluded" and has to reason backwards from two booleans to
work out whether that is a flag they can raise or a problem in their scoring file.

## Scope

Add a nullable `exclusion_reason` column, NULL unless `match_summary = 'excluded'`:

| `exclusion_reason` | meaning | set by | suppressed by |
|---|---|---|---|
| `multiallelic` | target carries more than one alternative allele | `match_variant`, before any allele comparison | `--keep_multiallelic` |
| `ambiguous` | a complemented reading of the score variant also matched | `match_variant`, per alternative allele | `--keep_ambiguous` |
| `collision` | another score variant of this accession resolved to the same target allele | the post-match pass, across rows | — |

`multiallelic` takes precedence over `ambiguous` where both apply. The existing code
order already gives this for free — `match_variant` tests `is_multiallelic and not
match_multiallelic` before matching any alleles (`lib/score/_matchvariants.py:190`).
`collision` is necessarily last, being the only reason that is not a property of a single
score variant.

Carry the column through `allele_match_table`, `score_log_table` and the summary's
`GROUP BY`, so each exclusion count names its own fix: raise `--keep_multiallelic`, raise
`--keep_ambiguous`, or go and look at the scoring file.

Two writers, deliberately. `match_variant` knows its own reason and sets it in the
`MatchResult` struct; the collision pass sets it in SQL, because a collision is not
visible from inside a single UDF call.

## Out of scope

Splitting `match_type_enum` into `EXCLUDED_AMBIGUOUS` / `EXCLUDED_MULTIALLELIC` /
`EXCLUDED_COLLISION`. `match_type` is derived from `MatchPriority.name`, and
`MatchPriority` is an ordering consumed by `max()` and by the ranking window
(`lib/score/_matchvariants.py:208`, lines 107-114). Three exclusion members need three
priority values whose relative order means nothing, encoding a false ranking into the
type that drives match selection. Rejected in ADR 0010.

Subcategorising `unmatched`. Nobody has asked for it and the reasons are not enumerable.

## Acceptance

- Every excluded row carries exactly one `exclusion_reason`; every non-excluded row
  carries NULL.
- A target that is both multiallelic and ambiguous, with both flags off, reports
  `multiallelic`.
- The summary log reports excluded counts split by reason, and `is_ambiguous` /
  `is_multiallelic` remain descriptive — a matched multiallelic variant still reports
  `is_multiallelic = true`.
- Raising `--keep_multiallelic` moves those variants out of the excluded counts and
  into matched or unmatched, without changing the other two reasons' counts.
