# 0030 — Only the first alternative allele can match, and ALT order decides the answer

**Refs:** ADR 0010.

## Why

`get_variant_match_priority` breaks after the first alternative allele when the variant
is multiallelic (`lib/score/_matchvariants.py:355-357`), but the ambiguous branch
`continue`s past that break (line 338). Two consequences, both bad.

**Most score variants at a multiallelic target cannot match at all.** Against
`REF=A ALTS=[C,G]`, on `main` today:

```
row1: effect A, other C  ->  REFALT    idx=0
row2: effect A, other G  ->  NO_MATCH  idx=None     <- should be REFALT idx 0
row3: effect C, other A  ->  ALTREF    idx=1
row4: effect G, other A  ->  NO_MATCH  idx=None     <- should be ALTREF idx 2
```

Only score variants pairing with `alts[0]` are reachable. `--keep_multiallelic` does not
help; it gates whether matching is attempted, not how far it gets.

**The answer depends on ALT column order.** `REF=A ALTS=[C,T]` against a score variant
`A/T` reports *unmatched* — alt `C` yields `NO_MATCH`, then the break fires. The same
target written `ALTS=[T,C]` reports *excluded* — alt `T` is ambiguous, which `continue`s
past the break, and `max()` then picks `EXCLUDED` over the later `NO_MATCH`. ALT order
carries no meaning, so neither should the result.

## Scope

Remove the break: match a score variant against every alternative allele and take the
best priority.

No tie-break policy is needed, and the reason is worth a comment rather than a
mechanism. At most one alternative allele can match a given score variant — with an
other allele, a forward match on `alt_i` and a complemented match on `alt_j` both require
`alt_i = complement(ref) = alt_j`; without one, the effect allele equals exactly one
distinct alternative allele. `max()` therefore never arbitrates between two real matches,
only between a match and the non-matches around it.

Unify ambiguity while the loop is open, because the two changes touch the same lines and
the pairwise check cannot survive the first one intact. A score variant is ambiguous when
**a complemented reading of it also matches an allele the target carries**:

| target | score variant | forward | complemented | verdict |
|---|---|---|---|---|
| `A / T` | `A/T` | idx 0 | idx 1 | ambiguous — palindrome, as before |
| `A / C` | `A/C` | idx 0 | — | resolvable |
| `A / C` | `C` (no other allele) | idx 1 | — | resolvable |
| `A / C,G` | `C` (no other allele) | idx 1 | idx 2 | ambiguous — only reachable multiallelically |
| `A / T,C` | `C/A` | idx 2 | — | resolvable despite the palindromic pair at the site |

This is a strict generalisation of `is_biallelic_variant_ambiguous`, so that function and
the `match_ambiguous` pre-check are deleted and ambiguity is computed from whether a
second reading matched.

Row four is the case that motivates it and the only place a strand flip changes which
dosage gets counted. `classify_match` attempts no complemented reading when the other
allele is NULL (lines 304-309), which is harmless biallelically — the complemented effect
allele is not at the site, so nothing can be confused with it. At `REF=A ALTS=[C,G]` both
readings land on real alleles and the forward one silently wins.

## Out of scope

Two score variants colliding on one target allele (ticket 0027). That is the opposite
direction — many score variants onto one allele, rather than many candidate alleles for
one score variant — and prioritisation cannot fix it.

## Acceptance

- All four rows of the table above match, with indices 0, 0, 1 and 2 respectively.
- `REF=A ALTS=[C,T]` and `REF=A ALTS=[T,C]` against `A/T` give the same result.
- `is_biallelic_variant_ambiguous` is gone, and the biallelic palindrome cases it
  covered still report `EXCLUDED` with `is_ambiguous = true`.
- A no-other-allele score variant is excluded when the target carries both the effect
  allele and its complement, and matched when it does not.
- **Expected log change, not a regression:** ambiguity is currently computed from
  `(ref, alt)` before anyone checks whether the score variant matched, so a score
  variant of `G/C` at a palindromic `A/T` target — matching nothing — is logged
  `EXCLUDED` today and becomes `NO_MATCH`. No score changes; `match_rate` sums only
  `matched` (`lib/score/_matchlog.py:137`). Update the fixtures that encode the old
  split rather than preserving it.
