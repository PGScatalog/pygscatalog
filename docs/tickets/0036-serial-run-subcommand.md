# 0036 — `pgsc_calc run`: execute the four stages in sequence

**Refs:** ADR 0013, ADR 0005, ADR 0006.

## Why

Nextflow needs a JVM and its own configuration, which is a real obstacle in locked-down
environments, and its fan-out is overhead rather than help when there is nothing to fan
out to. A **Serial run** (`CONTEXT.md`) closes that gap, and dogfoods the stage
interfaces this repo owns: today nothing here runs them end to end except a test fixture.

ADR 0006 rejected an in-process orchestrator. ADR 0013 explains why that rejection does
not reach a sequential runner, and what has to stay true for it not to become one.

## Scope

New subcommand `pgsc_calc run`, alongside the five that exist. Arguments are LOAD's
(`--target_genomes PATH:CHROM ...`, `--sampleset`, `--format`, `--score_paths`,
`--absence_semantics`, `--batch_size`, `--bgen_sample_file`), plus `--min_overlap`,
`--out_dir`, `--max_memory_gb` and `--threads`.

**LOAD runs once per target, into its own cache directory.** `GenotypeCache.shard` is
"the single leaf group this SCORE task holds" and `score` takes `--zarr_zip_file`
singular, so one zip per shard is a requirement, not a convention. The shard set is the
user's file list.

**Stages are invoked by composing argv and dispatching through the parser `main()`
builds** — `parser.parse_args(argv)` then `args.func(args)` — not by importing stage
functions and hand-building a `Namespace`. Argparse is then the drift check: a renamed
or newly-required stage argument fails naming the flag, instead of raising
`AttributeError` part-way through a stage that has already done an hour of work.

MATCH is invoked once per chromosome over every shard loaded under it; SCORE once per
shard; MATCH_COMBINE and SCORE_COMBINE once each. Layout:

```
out_dir/
  scores.txt.gz                                   <- copied from stages/combine/
  stages/
    load/0/genotypes.zarr.zip
    match/1/{matches.parquet, samples.parquet, scratch/}
    match_combine/{match_table.parquet, logs/, scratch/}
    score/0/{partial_scores.parquet, scratch/}
    combine/{scores.txt.gz, scratch/}
```

This tree is `run`'s output, not its scratch: ADR 0005 governs what a Nextflow
`publishDir` sees and still holds *within* each stage directory, which is what ticket
0037's tests keep asserting.

**Abort on the first stage failure.** SCORE_COMBINE sums the partial scores it is handed
and cannot know how many shards existed, so continuing past a failed SCORE returns
quietly low scores for every sample. Before invoking SCORE_COMBINE, assert one
`partial_scores.parquet` per shard — `run` is the only component that knows the expected
count. On failure, leave the tree in place and log the failing stage's argv, so that one
stage can be rerun by hand.

The stage sequence is a declared list that a test pins, so adding a stage is a choice
rather than an oversight (ADR 0013's residual risk).

## Out of scope

Resume. A failure costs the whole run; the recovery path is the logged argv above.

A dataset size guard. The limit is shard count, not sample or variant count, so any
threshold would measure the wrong quantity.

Per-stage resource limits. `--max_memory_gb` and `--threads` are declared once and apply
to every stage, where Nextflow sets them per task.

## Acceptance

- `pgsc_calc run` over the `grid_targets` shape reproduces the golden scores, matching
  what the same targets produce through the stages invoked individually.
- Every artefact in the tree above exists after a successful run, and `scores.txt.gz` is
  at the `out_dir` root.
- A SCORE failure aborts the run with a non-zero exit; SCORE_COMBINE is never invoked,
  and no `scores.txt.gz` is written.
- A missing partial score is refused by name before SCORE_COMBINE runs, not absorbed
  into a lower score.
- Removing a stage from the declared list fails a test.
