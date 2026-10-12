# 0001 — Nested cache hierarchy: /{sampleset}/{chromosome}/{filename}

**Blocks:** 0003, 0004, 0005, 0007
**Refs:** ADR 0003

## Why

The cache layout is hardcoded at depth 2, which is what actually prevents more than one file per chromosome — not the storage format. Adding a chromosome level makes the store self-describing and lets MATCH fan out per chromosome while SCORE fans out per shard.

## Scope

- `TargetGenome._zarr_group_path` → `f"{sampleset}/{chrom or 'ALL'}/{filename}"`. Unsplit whole-genome input uses the `ALL` sentinel so depth is uniform for every input shape.
- `get_variants_from_zarr` (`scorepipeline.py:352`) yields a path relative to the sampleset group instead of `group.basename`. The recursion already exists and is currently wasted.
- `unzip_zarr` (`load.py:238`) replaces the hardcoded `split("/")[0]` / `[0:2]` ancestor `zarr.json` extraction with a loop over all ancestors.
- `targetvariants`: replace the `filename` column with `zarr_group` holding the full relative path (`UKB/1/batch1.vcf.gz`). PK becomes `(variant_id, zarr_group)`. This deletes the `concat_ws('/', sampleset, filename)` reconstruction in `_weight_matrix.py:118` — the shard's identity is stored rather than rebuilt. Keep `sampleset` as a column for grouping.
- `write_dosage_regions` needs no change; zarr indexes `"a/b/c"` natively.

## Out of scope

Sample batching itself (0005), and any migration of alpha1 caches — they use the depth-2 layout and are discarded, not migrated. Re-running LOAD is expensive in a TRE but alpha2 is pre-release.

## Acceptance

- A sampleset with 2 chromosomes × 2 files per chromosome loads, zips, partially extracts one leaf, and scores.
- An unsplit whole-genome input lands under `ALL` and scores identically to before the change.
- Existing one-file-per-chromosome tests pass unmodified.
- No remaining code path assumes a fixed hierarchy depth (grep for `split("/")` and `basename` under `lib/`).
