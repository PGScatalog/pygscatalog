# 0015 — Delete the multi-shard dosage machinery

**Refs:** 0007, 0008, ADR 0003

## Why

SCORE now runs on one shard (0007), and `ScorePipeline.shard` raises if a task is staged with
more than one. `create_wide_weights_table` filters the match table to that single
`zarr_group`, so `wide_score_variants` holds exactly one distinct group, and both
`group_paths` and `df_groups` can only ever hold one entry.

The machinery for many is still there:

- `weight_mat_row_nr` is assigned by `row_number() OVER (ORDER BY zarr_group, target_row_nr,
  effect_allele_idx)` — a term in that ordering that now has one value
- `store_group_weight_arrays` loops over `group_paths`
- `store_results_in_zarr` has a `ContainsArrayError` branch appending to an existing
  `weight_matrix`, reachable only on a second group
- `write_dosage_regions` loops over `df_groups`, tracks a region offset, and asserts
  `start == min(weight_mat_row_nr)` to catch a misalignment that needs two groups to happen

The cost is not the dead code, it is that one invariant — *`weight_mat_row_nr` is the row
index of the dosage array* — is spread across five sites that must agree: the SQL window
function, the `group_paths` ordering, the `partition_by` dict ordering, the region slices, and
`calculate_score_statistics` using the row number to index the dosage array. Misalignment
produces a wrong score with no error and no NaN. Removing the many-group case makes the
invariant trivially true instead of checked.

## Scope

- Score one shard's weight matrix and dosage array end to end. Drop the group loops, the
  region offsets, the contiguity check and the append branch.
- `weight_mat_row_nr` becomes the position of a row in the shard's own ordering; drop
  `zarr_group` from the row-number ordering.
- Write the dosage and missing arrays whole rather than by region.

## Out of scope

Reintroducing multi-shard SCORE — 0007 and ADR 0003 put the fan-out in Nextflow. SCORE_COMBINE
is untouched: its reduction is valid because every `(sample, variant)` pair appears in exactly
one leaf (0008), which is a different invariant and stays.

## Acceptance

- Golden scores unchanged for every integration fixture: unsplit, chromosome-split, sample
  batches, and the chromosome × batch grid.
- No loop over zarr groups anywhere in the SCORE path.
- A shard where one scoring file matches REFALT and another ALTREF at the same position still
  produces two weight-matrix rows — the `effect_allele_idx` tie-break survives the change to
  the row ordering. This is the case the old ordering comment called out as "really
  important"; it needs a test that fails if the tie-break is dropped.
- The NaN sweep at the end of dosage calculation still runs and still fails a truncated write.
