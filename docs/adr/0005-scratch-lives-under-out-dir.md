# Scratch lives under `out_dir`; published artefacts stay at the root

## Status

Accepted

## Context

Before the alpha2 split, a stage's `out_dir` carried four meanings at once: the cache root (`genotypes.zarr`), the database location (`pgs.db`), the scratch store for the dosage array (`calculate_scores` opened a `LocalStore(out_dir)` and wrote `pgs/{sampleset}` plus a top-level `zarr.json` into it), and the directory the stage's published artefact was written to. A task's outputs were therefore not distinguishable from its scratch, and the dosage array — shaped `(n_score_variants, n_samples)`, the largest thing SCORE creates — landed in the published directory.

This is not purely an internal tidy-up. Each stage is its own Nextflow process ([ADR 0001](0001-nextflow-fanout-for-split-apply.md)), so `out_dir` is a task work directory and what leaves it is decided by the process definition in the separate `pgsc_calc` repo. Renaming or relocating a published artefact breaks wiring this repo does not control, which is what makes the layout a contract rather than a preference.

Two alternatives were considered.

A separate `--scratch_dir` argument would let a TRE put the dosage array on task-local storage instead of the published mount — which matters when that mount is network storage under quota, as it often is in a TRE. It was rejected for alpha2 because it costs a new argument on all four stages plus a `pgsc_calc` change, for a benefit nobody has measured yet.

Leaving the layout alone and documenting it is defensible, because the Nextflow processes name their outputs explicitly rather than globbing, so nothing is actually mis-published today. It was rejected because it leaves any future `publishDir` glob having to exclude `genotypes.zarr` and `pgs.db` by name, and leaves "what did this stage produce?" unanswerable from the directory alone.

## Decision

Scratch goes in a subdirectory of `out_dir`. Published artefacts stay at the `out_dir` root, with the same filenames.

```
out_dir/
  matches.parquet          <- declared output
  samples.parquet          <- declared output
  scratch/
    pgs.db
    genotypes.zarr/
    pgs/{sampleset}/dosage
```

Every artefact a Nextflow process declares — `matches.parquet`, `samples.parquet`, `match_table.parquet`, `partial_scores.parquet`, `scores.txt.gz`, `logs/` — keeps its name and its place. No new CLI argument.

## Consequences

- `pgsc_calc` needs no change: the artefacts it declares are where they were.
- A stage's published outputs are enumerable from the directory: its top level is exactly the declared artefacts plus `scratch`. This is asserted by a test rather than left as a convention.
- A future `publishDir` glob has one directory to ignore rather than a list of scratch filenames to name.
- **This decision separates naming, not storage.** The dosage array still sits on the same filesystem as the published output. If a TRE needs it elsewhere, `--scratch_dir` is the answer and remains strictly additive: defaulting it to `out_dir/scratch` keeps every current caller working.
- The scratch path is not a stable interface. Nothing outside a stage should read `out_dir/scratch`, and one integration test that reads the dosage array had to follow it there — a signal that any future test doing the same is reaching past a stage's interface.
