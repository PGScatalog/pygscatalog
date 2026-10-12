# `pgsc_calc run` executes the four stages sequentially for convenience

## Status

Accepted

## Context

[ADR 0006](0006-one-module-per-stage-no-base-class.md) rejected an orchestrator
that runs all four stages in one process. It gave two reasons.

The first reason concerned scaling.
[ADR 0001](0001-nextflow-fanout-for-split-apply.md) assigns fan-out to Nextflow.
A dask-distributed scheduler needs connections between nodes, which TREs often
forbid. The second reason concerned maintenance. Another orchestrator would define
the stage graph again and could diverge from `main.nf`.

The first reason does not apply to a sequential runner. Sequential execution on
one machine adds no cluster, scheduler or network requirement. The maintenance
risk remains. It increases when the runner becomes a shipped feature instead of a
test fixture.

Demand changed. Nextflow requires a JVM and its own configuration, which can block
use in restricted environments. Fan-out adds overhead when no other execution
resources exist.

There is also a reason to test the entry point that users receive. Currently,
only `run_pipeline` in `tests/integration/conftest.py` exercises the repository's
stages end to end. This helper already defines the graph a second time. ADR 0006
accepted that duplication because the helper did not ship.

## Decision

A fifth subcommand, `pgsc_calc run`, loads each target genome. It then executes
MATCH → MATCH_COMBINE → SCORE → SCORE_COMBINE sequentially in one process.
`CONTEXT.md` calls this a **Serial run**. It provides convenience, not another
scaling mechanism. `pgsc_calc` remains the supported path. The serial entry point
is expected to serve a small minority of users.

**The change deletes `run_pipeline`, so it does not add another stage graph.**
Integration tests drive `pgsc_calc run`. The repository still has one stage-graph
definition, and tests exercise the shipped definition.

**The runner composes argv and dispatches through the parser that `main()`
builds.** It does not import stage functions and construct a `Namespace` manually.
Argparse checks that `run` uses each stage's documented arguments. A renamed or
newly required argument produces an error that names the flag. It does not cause
an `AttributeError` during the stage.

**Stages exchange the same disk artefacts as they do under Nextflow.** Passing
polars frames between stages in memory was rejected. The serial runner exercises
the artefact contract that binds this repository to `main.nf`. A serial run
also archives and extracts the cache again, as a Nextflow task does. This keeps
the LOAD → MATCH/SCORE boundary under test.

Each stage has a directory under `out_dir/stages/`. The runner copies the final
`scores.txt.gz` to the `out_dir` root:

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

This tree is `run`'s output, not its scratch.
[ADR 0005](0005-scratch-lives-under-out-dir.md) governs what a Nextflow
`publishDir` sees. It still applies within each stage directory. It does not
govern a serial run's root, because a serial run is not a Nextflow task.

**`run` aborts at the first stage failure.** SCORE_COMBINE sums the supplied
partial scores. It cannot know how many shards existed. Continuing after a failed
SCORE would silently reduce every sample's score. The pipeline refuses that
failure elsewhere too.

Before SCORE_COMBINE, `run` asserts that it has one partial score per shard. Only
`run` knows the expected count.

**Resume is outside scope.** A failure requires the whole run again. Outputs from
completed stages remain on disk. The runner logs the failing stage's argv. A user
can run that stage manually against those outputs.

## Consequences

- A serial run's peak memory is the maximum across stages, not their sum. SCORE
  still holds one shard ([ADR 0003](0003-shards-are-zarr-leaf-groups.md)), and MATCH
  still sees one chromosome. Sequential execution increases elapsed time.
  **Shard count therefore limits serial use, rather than dataset size.** A few
  large shards suit it better than many small shards.

  There is deliberately no size guard. A size threshold would measure the wrong
  quantity.
- Peak memory equals the largest stage only if each stage releases its resources.
  Leaked stores or open connections can make usage accumulate across stages.
  Nextflow avoids this accumulation by ending each process. A test comparing peak
  RSS with one SCORE task would be useful.
- The runner declares `--max_memory_gb` and `--threads` once and applies them to
  every stage. Nextflow sets these resources per task.
- **Nothing detects a fifth stage.** Argparse catches changed arguments. It cannot
  detect a stage added to `main.nf` but omitted from `run`. Every test here would
  continue to pass. This is an accepted risk.

  The only mitigation is a declared stage list that `run` iterates and a test
  checks. Adding a stage then requires an explicit choice. A coverage gate was
  rejected. The new stage's own tests would cover its lines while the serial
  runner still omitted it.
- Tests continue to assert the stage output layout. Only their prefix changes to
  `out_dir/stages/`. They do not use scratch paths, which ADR 0005 declares
  unstable.
