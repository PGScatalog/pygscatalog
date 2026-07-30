# 0012 — Give the Cache a module

**Blocks:** 0013, 0014, 0017
**Refs:** ADR 0003, ADR 0004

## Why

`CONTEXT.md` calls the Cache "the portable interface between stages", and alpha2 leans on it
harder than anything else: shards are defined over it, stages fan out at different levels of
it, and every invariant is checked against it. It is the only one of those concepts with no
module. Nine callers each carry their own copy of the hierarchy — `get_shard_attrs`,
`get_variants_from_zarr`, `check_cache_shape`, `ScorePipeline.sampleset` /
`sampleset_group` / `shard` / `samplesets` / `genotypes_root`, plus the traversals in 0013
and the path construction in 0014.

Two concrete costs:

- `check_cache_shape` recovers the chromosome with `path.split("/")[1]`. ADR 0003's stated
  consequence was that three of four code sites *delete* a hardcoded depth; this one grew a
  fresh one back, in the function the ADR introduced.
- `check_cache_shape`, `get_shard_attrs` and `get_variants_from_zarr` have no unit tests.
  The three leaf invariants (`check_absence_semantics`, `check_samples_disjoint`,
  `check_shards_are_rectangular`) are unit-tested against hand-built dicts, but the function
  that finds the shards and calls them is reachable only through the integration suite,
  because its interface takes a `zarr.Group`.

## Scope

- A module in `lib/cache` owning the hierarchy. Its interface answers what a caller actually
  wants to know: which shards a staged cache holds (path, samples, declared absence
  semantics), which chromosome it holds, and which single shard a SCORE task holds.
- The chromosome comes from a shard's position in the tree, not from slicing its path string.
- `check_cache_shape` and `get_shard_attrs` move behind it, along with `ScorePipeline`'s
  `sampleset`, `sampleset_group`, `shard`, `samplesets` and `genotypes_root`. The two
  different single-sampleset assertions (`NotImplementedError` in `sampleset` and a bare one
  in `samplesets`) collapse into one.

## Out of scope

The remaining reader traversals (0013) and the writer-side leaf path (0014) — this ticket
lands the seam, those two fold callers in behind it. Splitting `ScorePipeline` (0017).

Nothing about the hierarchy itself changes, and no invariant changes what it checks. This is
a move, not a redesign: if scores or error messages change, the ticket has gone wrong.

## Acceptance

- `check_cache_shape` is unit-testable from a fake set of shards, with no zarr store on disk,
  and gets a test for the case it currently has none for: a cache holding two chromosomes.
- No depth assumption about the hierarchy survives outside the module (grep `split("/")` and
  `basename` under `lib/`).
- One traversal of the shard tree, not four.
- Integration tests pass unmodified, including the golden scores.
