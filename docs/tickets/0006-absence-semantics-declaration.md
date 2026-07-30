# 0006 — Declare absence semantics per sampleset

**Refs:** ADR 0004

## Why

A variant absent from a target genome file means either "never assessed" (`unobserved`, the historical imputed-array assumption `min_overlap` is built on) or "assessed, and every sample here is homozygous reference" (`homref`, sparse-called WGS). The two produce byte-identical files, so this cannot be detected — only declared.

Declaring it now, while only one mode is implemented, is what keeps `homref` additive later. Retrofitting the declaration after the fact would mean changing the default behaviour of existing samplesheets.

## Scope

- Samplesheet column, one value per sampleset, defaulting to `unobserved`. Not per shard: mixed semantics within a sampleset means the shards are not slices of one callset, which makes them separate samplesets.
- Plumb through to MATCH so 0005's ragged check can consult it.
- Only `unobserved` is accepted in alpha2. `homref` is rejected with "not yet implemented", not silently ignored.
- Document that bgen is dense, so every bgen shard records every variant it assessed and the declaration cannot change bgen behaviour.

## Out of scope

The `homref` code path entirely: constant fill rows, union-derived allele indices, and the gate treating an assessed-elsewhere variant as matched. All strictly additive — see ADR 0004 for the design.

## Acceptance

- A samplesheet without the column loads and behaves exactly as today.
- `homref` produces a clear not-implemented error, not a wrong score.
- The value reaches MATCH and drives 0005's ragged check.

## Not in this repo

The samplesheet column itself lives in `pgsc_calc`'s `assets/schema_input.json`, which
sets `additionalProperties: false`, so the column has to be declared there before a
samplesheet can carry it. This repo exposes it as `pgsc_calc load --absence_semantics`
and records the declaration on each **leaf** zarr group, not the sampleset group: LOAD
is fanned out per shard and merging the per-shard zips is last-write-wins on an
ancestor `zarr.json`, so a sampleset-level attribute would silently swallow a
contradictory declaration instead of letting MATCH detect it.
