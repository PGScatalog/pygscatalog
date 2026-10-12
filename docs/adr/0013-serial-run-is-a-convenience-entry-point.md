# `pgsc_calc run` executes the four stages in sequence, as a convenience entry point

## Status

Accepted

## Context

[ADR 0006](0006-one-module-per-stage-no-base-class.md) rejected "an orchestrator that runs all four stages in one process", on two grounds: it contradicts [ADR 0001](0001-nextflow-fanout-for-split-apply.md), where fan-out belongs to Nextflow because a dask-distributed scheduler needs inter-node reachability that TREs often forbid; and it would be a second definition of the stage graph, free to diverge from `main.nf`.

The first ground does not reach a *sequential* runner. It is an argument against a second scaling mechanism, and running stages one after another on one machine introduces no cluster, no scheduler, and no network requirement. The second ground stands, and gets sharper the moment such a runner is shipped rather than confined to a test fixture.

What has changed is the demand. Nextflow needs a JVM and its own configuration, which is a real obstacle in locked-down environments, and its fan-out is overhead rather than help when there is nothing to fan out to. There is also a dogfooding argument: this repo's stages are currently exercised end-to-end only by `run_pipeline` in `tests/integration/conftest.py`, which is itself the second definition ADR 0006 tolerated on the grounds that it was not shipped.

## Decision

A fifth subcommand, `pgsc_calc run`, executes MATCH → MATCH_COMBINE → SCORE → SCORE_COMBINE in sequence in one process, after LOADing each target genome. `CONTEXT.md` names this a **Serial run**. It is a convenience entry point, not a scaling mechanism: `pgsc_calc` remains the supported path, and this is expected to serve a small minority of users.

**It is not a second stage graph, because `run_pipeline` is deleted.** The integration tests drive `pgsc_calc run`, so the count of stage-graph definitions in this repo stays at one, and the one that ships is the one the tests exercise.

**Stages are invoked by composing argv and dispatching through the parser `main()` builds**, not by importing stage functions and hand-building a `Namespace`. Argparse then validates that `run` still speaks each stage's documented argument language: a renamed or newly-required argument fails with a clean error naming the flag, rather than an `AttributeError` part-way through a stage.

**Stages exchange the same on-disk artefacts they exchange under Nextflow.** Passing polars frames between stages in memory was rejected: the artefact contract is what binds this repo to `main.nf`, and it is the thing worth dogfooding. A serial run also archives and re-extracts the cache like a Nextflow task does, keeping the LOAD → MATCH/SCORE seam covered.

Each stage gets its own directory under `out_dir/stages/`, and the final `scores.txt.gz` is copied to the `out_dir` root:

```
out_dir/
  scores.txt.gz
  stages/
    load/0/genotypes.zarr.zip
    match/1/{matches.parquet, samples.parquet, scratch/}
    match_combine/{match_table.parquet, logs/, scratch/}
    score/0/{partial_scores.parquet, scratch/}
    combine/{scores.txt.gz, scratch/}
```

This tree is `run`'s output, not its scratch. [ADR 0005](0005-scratch-lives-under-out-dir.md) governs what a Nextflow `publishDir` sees and continues to hold *within* each stage directory; it does not govern the root of a serial run, which is not a Nextflow task.

**`run` aborts on the first stage failure.** SCORE_COMBINE sums the partial scores it is handed and cannot know how many shards existed, so continuing past a failed SCORE would produce quietly low scores for every sample — the failure mode this pipeline refuses elsewhere. `run` asserts it has one partial score per shard before invoking SCORE_COMBINE, because it is the only component that knows the expected count.

**Resume is out of scope.** A failure costs the whole run. The completed stages' outputs are left on disk and the failing stage's argv is logged, so that stage can be rerun by hand against them.

## Consequences

- A serial run's peak memory is the maximum over single stages, not their sum, because SCORE still holds one shard ([ADR 0003](0003-shards-are-zarr-leaf-groups.md)) and MATCH still sees one chromosome. What a serial run gives up is wall-clock. **Its limit is therefore shard count, not dataset size** — big-but-few shards are a good fit, and many-small-shards is the miserable case. There is deliberately no size guard, because any threshold would be measuring the wrong quantity.
- One process means "peak = max over stages" only if each stage releases what it held. A leaked store or an open connection turns max into sum, which Nextflow gets right for free by ending the process. Worth a test that pins peak RSS against a single SCORE task.
- `--max_memory_gb` and `--threads` are declared once and apply to every stage, where Nextflow sets them per task.
- **Nothing detects a fifth stage.** Argparse catches a changed argument, but a stage added to `main.nf` and never wired into `run` will keep passing every test here. This is an accepted risk, mitigated only by `run` iterating a declared list of stages that a test pins, so adding one is a choice rather than an oversight. A coverage gate was considered and rejected: a new stage's own tests would cover its lines, so the gate would report healthy for precisely this failure.
- The tests that assert stage output layout keep asserting it — only their prefix moves, to `out_dir/stages/`. They are not re-pointed at a scratch path, which ADR 0005 declares unstable.
