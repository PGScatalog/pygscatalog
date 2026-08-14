# 0040 — Sex check: observed hemizygous fraction against the declaration

**Refs:** ADR 0014, ADR 0007, ADR 0009, ADR 0003.

**Depends on:** 0039.

## Why

Nothing in the pipeline notices that the sample labelled female is hemizygous across
chrX. A sample swap, a mislabelled cohort, or a sex map keyed on the wrong ID space all
score cleanly and silently — and since ADR 0007, X variants contribute to those scores.

ADR 0007 stays intact: ploidy is still observed, the declaration is still never used to
derive it. This ticket only compares the two and stops when they disagree.

## Scope

Runs inside LOAD, per shard. **No new stage**: `CONTEXT.md` guarantees every
`(sample, variant)` pair appears in at most one shard, and multiple files under one
chromosome are sample batches with disjoint samples (ADR 0003), so one chrX shard holds
every cached X call for its samples. There is nothing to combine.

**Informative variants.** A cached X variant is informative when at least one call in the
shard is hemizygous. A variant diploid in every sample is excluded: it is either
pseudoautosomal or a cohort with no hemizygous samples in it, and nothing tells those
apart without the build-specific PAR coordinates ADR 0007 declined to ship.

This exclusion is not a refinement. A scoring file whose X variants are mostly
pseudoautosomal makes a man diploid at most of them, landing his fraction near 0 and
reading him as female — without the rule the check is wrong, not weak.

Ploidy comes from the array ADR 0007 already writes — `(gt != NO_ALLELE).sum(axis=2)` —
so this adds no format-specific code and works identically on both readers.

**The fraction.** Per sample, over informative variants, the proportion of its *called*
genotypes that are hemizygous. Classify on the **majority**: declared male expects 1,
declared female expects 0. No allele frequencies, no F coefficient, no LD pruning — ADR
0014 records the four reasons F does not survive contact with this cache.

**Below ten informative variants the check does not conclude**, and warns loudly naming
the count. Hemizygosity is deterministic rather than estimated, so the floor is tolerance
to bad calls, not statistical power: at one informative variant a single miswritten
no-call inverts the fraction and stops a run over a sample that is fine. Not concluding
is a warning rather than a failure — refusing every cohort whose PGS lacks X coverage
would make 0039's mandate self-defeating.

**No chrX among the targets: skip and warn.** The map was still mandatory and goes unread.

**`n_informative == 0` yields a NULL fraction, never `0.0`.** Computed naively, `0/0` is
`0.0`, which reads as diploid, which reads as female-consistent: every woman passes
silently and every man is flagged, on a check that never ran. Same class of error as
`sentinel_to_nan` in ADR 0007, and it needs a test that fails loudly rather than a
comment.

**A mismatch raises. `0` is the override.** Where a sample is declared male or female and
the majority of its informative calls says otherwise, LOAD raises **naming every
mismatching sample at once** — failing on the first would make a cohort with four hundred
mismatches four hundred round trips. Users exempt individual samples by declaring them
`0`, which is the documented recovery path, not a workaround.

**Raise after `zip_zarr`, not at the point of discovery.** The sex map influences nothing
that gets cached, so a cache written before the exception is valid and complete: the user
fixes the map, re-runs with `--zarr_zip_file`, and `get_unique_positions` returns nothing
(`cli/load.py:245`), so the retry costs a check rather than a re-load. This looks
backwards — LOAD deliberately finishes work it is about to reject — so it needs a comment
citing ADR 0014, or someone will tidy it into a fast failure and make every override
iteration re-read chrX.

**Reporting.** A per-sample TSV (`sample_id`, `declared`, `hemizygous_fraction`,
`n_informative`, `outcome`), written whether or not the run raises. Two outcomes for the
MVP: consistent and mismatch.

**Every run logs how many samples were exempted by declaration alongside how many were
checked**, unconditionally. This is the mitigation for the one failure mode ADR 0014
calls worse than the status quo: a user facing four hundred hard failures sets the whole
column to `0` and re-runs, and the pipeline exits cleanly having verified nothing. A run
that checked 3 of 500,000 samples has to say so in the same breath as it says it passed.

## Out of scope

**chrY as a second signal.** Only bgen can express the nullploidy that would make a
woman's chrY informative (ADR 0007), so it would make the check's strength
format-dependent. X carries it.

**A third outcome for intermediate fractions** (mosaicism, or a broken file). Collapses
into mismatch for the MVP; the fraction is in the TSV so it stays diagnosable. Worth
splitting later — it is not a finding about sex.

**A warning threshold on the exemption rate.** The unconditional count above is enough
until someone actually blanket-zeroes a column.

Distinguishing a genuinely unknown sex from a deliberate exemption. Both are `0`.

## Acceptance

- A mixed-sex chrX shard whose declarations agree with the calls loads and scores exactly
  as it does today.
- One sample's declaration flipped in the map fails LOAD, naming that sample.
- Two flipped samples are both named by a single failure.
- Flipping the same sample to `0` loads cleanly, and the log reports one exempted sample.
- A rerun after a failure, with a corrected map and `--zarr_zip_file`, re-reads no
  genotypes.
- A shard whose informative X variants number fewer than ten warns, does not conclude,
  and does not fail — including when a declaration is wrong.
- A shard whose cached X variants are all uniformly diploid reports zero informative
  variants and a NULL fraction, and passes no sample as female.
- An all-female shard cannot be checked and says so, rather than reporting every sample
  consistent.
- The TSV is written on both the passing and the failing path.
