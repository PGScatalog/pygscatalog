# 0014 — Build the leaf path once

**Depends on:** 0012
**Refs:** ADR 0003, 0002

## Why

Two functions independently build the same rule:

- `TargetGenome._zarr_group_path` → `{sampleset}/{chrom or 'ALL'}/{filename}`, which is where
  a shard actually gets written
- `AnnotatedGenomePath.leaf_path` → `{chrom or 'ALL'}/{filename}`, which is what 0002's
  uniqueness validation checks

So the validation guards a string the writer does not use. They agree today by coincidence of
being written together. If they drift, the guard passes and the failure it exists to prevent
comes back — and ADR 0003 named that failure precisely: two shards of the same chromosome
sharing a basename collide when the per-shard archives are merged, one leaf silently
overwriting the other. No error, no NaN, just a sampleset quietly missing genotypes.

## Scope

- One function derives a shard's leaf path. Both the writer and 0002's validation call it.
- Validation stays pure path arithmetic with no file access, as 0002 requires.

## Out of scope

Changing the rule itself, including the `ALL` sentinel and the decision not to derive the
group name from the parent directory — 0002 records why that does not work (Nextflow stages
`path` inputs flat, `/` is a zarr path separator).

## Acceptance

- A test asserting that, for the same input, the path 0002 validates and the path
  `TargetGenome` writes to are the same string.
- 0002's acceptance cases still hold: two inputs resolving to one leaf fail before LOAD and
  name both; the same basename under different chromosomes is accepted; two whole-genome VCFs
  with disjoint samples land under `ALL`.
