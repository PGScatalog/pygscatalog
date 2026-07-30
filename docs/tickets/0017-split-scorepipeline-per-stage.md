# 0017 — Split `ScorePipeline` into four stage modules

**Depends on:** 0012, 0016
**Refs:** ADR 0001, ADR 0005, ADR 0006

**Decision made** — recorded as [ADR 0005](../adr/0005-scratch-lives-under-out-dir.md) (scratch layout) and [ADR 0006](../adr/0006-one-module-per-stage-no-base-class.md) (the shape of the split).

## Why

`ScorePipeline` is one module whose interface is the union of four stages' needs. Only
`load_scores` is used by more than two stages. Three parts of the interface are
stage-conditional, and none of that is visible in a type signature:

- `shard` raises for three of the four stages
- `sampleset` reads the cache for MATCH and SCORE, but the two combine stages have no cache to
  read and must pass it to the constructor instead
- `minimum_samples_for_impute` is only meaningful to SCORE

`out_dir` carries four meanings at once: the cache root (`genotypes.zarr`), the database
location (`pgs.db`), the scratch store for the dosage array (`calculate_scores` opens a
`LocalStore(out_dir)` and writes `pgs/{sampleset}` plus a top-level `zarr.json` into it), and
the directory a stage's published artefact is written to. A task's outputs are therefore not
distinguishable from its scratch, and the dosage array — the largest thing SCORE creates —
lands in the published directory.

## Decided

Scratch goes in a subdirectory of `out_dir`; published artefacts stay at the `out_dir` root,
with the same filenames. No new CLI argument.

```
out_dir/
  matches.parquet          <- declared output, stays at root
  samples.parquet          <- declared output, stays at root
  scratch/
    pgs.db
    genotypes.zarr/
    pgs/{sampleset}/dosage
```

Every artefact a Nextflow process declares — `matches.parquet`, `samples.parquet`,
`match_table.parquet`, `partial_scores.parquet`, `scores.txt.gz`, `logs/` — keeps its name and
its place, so the wiring in `pgsc_calc` doesn't move. `publishDir` globs that would otherwise
have to exclude `genotypes.zarr` and `pgs.db` by name now only have to ignore one directory.

## Scope

- Four modules, one per stage, each taking only what it needs, over the cache seam from 0012
  and the database seam from 0016.
- Stage-conditional parts become constructor arguments of the stage that has them, so no
  member raises depending on which stage constructed it.
- Scratch and published artefacts separated per the decision above.

## Out of scope

- Artefact filenames and schemas. They are a real weakness — `matches.parquet` and
  `match_table.parquet` are structurally interchangeable, so feeding SCORE the ungated one
  scores accessions that failed `min_overlap` with no error — but that is its own ticket and
  changing both at once makes the diff unreadable.
- The stage graph, which stays in Nextflow (ADR 0001). This ticket reshapes the Python
  entrypoints each task calls, not who calls them.

## Acceptance

- No member of any stage module raises because of which stage constructed it.
- Each stage's CLI constructs only its own module.
- A stage's published outputs are enumerable without listing its scratch.
- Integration tests pass unmodified, and the four-stage run still produces the golden scores.
