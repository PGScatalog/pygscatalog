# Store scratch data under `out_dir` and keep published artefacts at its root

## Status

Accepted

## Context

Before the alpha2 split, each stage used `out_dir` for four purposes:

- The cache root, `genotypes.zarr`.
- The database, `pgs.db`.
- The scratch store for dosage arrays.
- The stage's published artefacts.

`calculate_scores` opened `LocalStore(out_dir)` and wrote `pgs/{sampleset}` plus a top-level `zarr.json` there.
The directory did not distinguish outputs from scratch data.
The dosage array had shape `(n_score_variants, n_samples)` and was the largest object that SCORE created.
It occupied the published directory.

This layout affects more than internal organisation.
Each stage runs as a separate Nextflow process, as described in [ADR 0001](0001-nextflow-fanout-for-split-apply.md).
Its `out_dir` is therefore a task work directory.
The process definition in the separate `pgsc_calc` repository determines which files leave that directory.

Renaming or moving a published artefact breaks connections outside this repository's control.
The layout is therefore a contract, not a preference.
We considered two alternatives.

A separate `--scratch_dir` argument would let a TRE place the dosage array on task-local storage instead of the published mount.
This matters when the mount uses network storage with a quota, as TRE mounts often do.
We rejected it for alpha2 because it requires a new argument on all four stages and a change in `pgsc_calc`.
Nobody measured the benefit before this decision.

Retaining and documenting the existing layout was also defensible.
Nextflow processes named outputs explicitly instead of using glob patterns, so they published no incorrect files at the time.
However, a future `publishDir` glob would need explicit exclusions for `genotypes.zarr` and `pgs.db`.
The directory alone would also fail to show what the stage produced.
We rejected this alternative for those reasons.

## Decision

Scratch data goes into a subdirectory of `out_dir`.
Published artefacts retain their filenames at the `out_dir` root.

```
out_dir/
  matches.parquet          <- declared output
  samples.parquet          <- declared output
  scratch/
    pgs.db
    genotypes.zarr/
    pgs/{sampleset}/dosage
```

Every artefact declared by a Nextflow process retains its name and location:

- `matches.parquet`
- `samples.parquet`
- `match_table.parquet`
- `partial_scores.parquet`
- `scores.txt.gz`
- `logs/`

The CLI gains no argument.

## Consequences

- `pgsc_calc` needs no change because its declared artefacts retain their locations.
- A stage's top level contains exactly its declared artefacts and `scratch`.
  A test checks this layout instead of relying on convention.
- A future `publishDir` glob needs to exclude one directory rather than several scratch filenames.
- **This decision separates names, not storage.**
  The dosage array remains on the same filesystem as the published output.
  A TRE that needs another location would require `--scratch_dir`.
  That addition remains compatible with existing callers if its default is `out_dir/scratch`.
- The scratch path is not a stable interface.
  Nothing outside a stage should read `out_dir/scratch`.
  One integration test read the dosage array and needed the new location.
  Any future test that reads scratch data likewise uses information outside the stage's interface.
