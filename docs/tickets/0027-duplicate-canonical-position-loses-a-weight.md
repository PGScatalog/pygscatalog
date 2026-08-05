# 0027 — Two scoring variants at one canonical position silently lose a weight

**Refs:** ADR 0008. Found while implementing 0019, whose last acceptance criterion asked
only that this be recorded.

> **Superseded in part by [ADR 0010](../adr/0010-multiallelic-matched-per-alt-bgen-stays-biallelic.md).**
> The acceptance criteria below stand. The **Scope** section's proposed
> `(accession, chr_name, chr_position)` uniqueness check on `score_variant_table` does
> not: once multiallelic variants are matched per alternative allele, an accession
> legitimately holds two rows at one position (`A/C` and `A/G` match different
> alternative alleles), and that key would reject valid scoring files. The check belongs
> after matching, keyed on `(accession, target_row_nr, effect_allele_idx)`. ADR 0010 also
> narrows the inflated-match-rate wrong: two score variants matching two different
> alleles of one variant are two variants found, not one locus counted twice.

## Why

0019 collapses `XY` into `X`, so a scoring file carrying an `XY` variant and an `X`
variant at the same position now has two rows at one canonical position. Run through the
real stages, that produces neither a collapse nor an error:

1. `score_variant_table` holds **two rows, both `X`, same position**, with their two
   different weights. `is_duplicated` stays `False` on both — pgscatalog-format computed
   it before `XY` became `X`, and nothing recomputes it afterwards.
2. `allele_match_table` matches **both** rows against the same `target_row_nr`. The PK is
   `(sampleset, accession, row_nr, zarr_group)`, so nothing dedups them: one target locus
   counts as two matched scoring variants in the `min_overlap` numerator.
3. `create_wide_weights_table`'s PIVOT groups on
   `(sampleset, zarr_group, is_recessive, is_dominant, target_row_nr, effect_allele_idx)`
   — identical for both rows — and aggregates with `coalesce(any_value(effect_weight), 0)`.
   The result is **a single weight-matrix row carrying one of the two weights**.

So the `GROUP BY` prevents a double count in the score itself, but **which weight
survives is unspecified** — both rows tie on the insert's `ORDER BY m.target_row_nr`.
A weight is discarded silently.

This is not confined to `XY`. Any two scoring variants that canonicalise to the same
position hit the same path; `XY`/`X` is just the first way to reach it.

## Scope

Decide what a duplicate canonical position means, and make the pipeline say so instead
of picking one arbitrarily. `score_variant_table`'s `PRIMARY KEY(row_nr, accession)` is
the natural place to add a `(accession, chr_name, chr_position)` uniqueness check if the
answer is "reject".

Note there are two separate wrongs here and a fix should address both: the discarded
weight, and the inflated match-rate numerator that counts one locus twice.

## Out of scope

Recomputing `is_duplicated` inside pgscatalog-format. The flag is computed upstream of
canonicalisation; making it correct here does not require making it correct there.

## Acceptance

- Two scoring variants at the same canonical position no longer resolve to an
  unspecified choice between their weights.
- The match rate counts one target locus once.
- A scoring file with a genuine `XY` and `X` pair at one position produces the same
  outcome on every run.
