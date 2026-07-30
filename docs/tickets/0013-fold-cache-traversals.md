# 0013 — Fold the remaining cache traversals behind the Cache module

**Depends on:** 0012
**Refs:** ADR 0003

## Why

Four functions walk the same tree with four different ideas of what a leaf is:

- `get_sample_ids` walks the subtree looking for a `samples` attribute and requires every
  shard to agree
- `get_n_samples` calls it only to take a length
- `get_variants_from_zarr` recurses separately, looking for a child group named `meta`
- `get_position_df` opens `{group.path}/meta` directly, as does the writer in
  `targetvariants`

Each is correct today. The cost is that "what a shard looks like in the cache" is asserted in
four places, so a change to the hierarchy has four chances to be applied inconsistently — and
inconsistently is the dangerous way, because a traversal that silently finds no shards reads
as an empty cache rather than an error.

## Scope

- `get_sample_ids`, `get_n_samples`, `get_variants_from_zarr` and `get_position_df` become
  callers of 0012's module rather than traversers in their own right.
- The `meta` group path is constructed in one place, and both the reader and the writer in
  `targetvariants` use it.

## Out of scope

`get_sample_ids`' `NotImplementedError` on disagreeing sample lists. ADR 0003 settled this:
it stops being reachable without being removed, because a SCORE task is staged with one
leaf's subtree and there is nothing left to disagree with. Move the behaviour, keep it, don't
argue with it.

## Acceptance

- One traversal of the shard tree in the codebase.
- `get_sample_ids` behaves identically for a cache with one leaf and for a cache whose leaves
  disagree — same return value, same exception.
- The `meta` path appears once.
- Integration tests and golden scores unchanged.
