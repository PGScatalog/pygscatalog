# 0003 — MATCH stage, fanned out per chromosome

**Depends on:** 0001
**Blocks:** 0004
**Refs:** ADR 0003

## Why

Matching currently runs inline at the start of a single `pgscatalog-calc score` invocation over the whole unified store. Splitting it out is what lets SCORE be gated on accessions that passed `min_overlap`, so rejected accessions never cost any scoring compute.

Fanning out per *chromosome* rather than per shard is the load-bearing choice: a chromosome-level task holds every shard's `meta` for that chromosome, so it matches against a complete variant list and no shard has a partial view. It reads 1D variant metadata only, never the 3D genotype arrays, so it does not scale with sample count.

## Scope

- New CLI: `pgscatalog-calc match`, taking one chromosome's leaves and the scoring files.
- Reuse `_load_variants_from_zarr` and `update_match_table` over the chromosome's leaves rather than the whole sampleset.
- Emit a partial match table carrying the shard's `zarr_group`, plus the chromosome's assessed variant set (currently the union of the leaves' variants).
- Per-shard `effect_allele_idx` resolution is unchanged — allele indices are already shard-local (`_weight_matrix.get_variant_metadata` keyed on `group_path`), which is what makes differing ALT orderings between shards safe.

## Out of scope

The `min_overlap` decision, which needs every chromosome (0004). MATCH emits counts; it does not gate.

## Acceptance

- Matching one chromosome's shards produces the same match table rows as today's whole-store run, restricted to that chromosome.
- Peak memory is flat as the number of samples per shard grows (no genotype array reads).
- Two shards under one chromosome with differing ALT orderings at a shared position both resolve their own `effect_allele_idx` correctly.
