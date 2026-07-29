# 0008 — SCORE_COMBINE: reduce partial score tables

**Depends on:** 0007
**Refs:** ADR 0003

## Why

Once SCORE is fanned out, no task produces final output. Combining is a single duckdb reduction.

## Scope

- New CLI: `pgscatalog-calc combine`. Union every shard's partial `score_table`, `GROUP BY (sampleset, accession, sample_id)`, `SUM` the numeric columns.
- No axis-awareness anywhere in this stage. The reduction is valid because **every `(sample, variant)` pair appears in exactly one shard** — variant-axis shards contribute disjoint variants for a shared sample, sample-axis shards contribute disjoint samples. Both reduce identically. (ADR 0002 justified this as "correct because every shard shares the same samples", which is narrower than what actually holds, and is why this stage needs no change to support batching.)
- Existing export paths (`export_scores`, match log exports) hang off the combined table.

## Acceptance

- Chromosome-split input reproduces today's scores.
- Sample-batched input reproduces the same scores as the same cohort in one file.
- Mixed chromosome × batch grid reproduces both.
- `n_matched` per sample equals the genome-wide count, not a per-shard count.
