# One module per stage; no stage base class, and no in-process orchestrator

## Status

Accepted

## Context

Scoring's four stages ([ADR 0001](0001-nextflow-fanout-for-split-apply.md), [ADR 0003](0003-shards-are-zarr-leaf-groups.md)) were driven by a single `ScorePipeline` class whose interface was the union of all four stages' needs. Only `load_scores` was used by more than two stages. Three parts of the interface were stage-conditional, and none of that was visible in a type signature:

- `shard` raised for three of the four stages
- `sampleset` was read from the cache by MATCH and SCORE, but had to be passed to the constructor by MATCH_COMBINE and SCORE_COMBINE, which have no cache to read
- `minimum_samples_for_impute` meant something only to SCORE

So a caller had to know which stage it was in order to know which members were legal, and the class grew a check for each of those facts instead of the facts being expressed by construction.

Two shared seams had meanwhile been separated out: the Cache (`GenotypeCache`, one module owning the zarr hierarchy) and the score database (`ScoreDatabase`, one module owning connections and the task's resource limits).

## Decision

Four concrete modules, one per stage, named after the stages in `CONTEXT.md`: `MatchStage`, `MatchCombineStage`, `ScoreStage`, `ScoreCombineStage`. Each takes only what it needs, over the two shared seams. Stage-conditional parts become constructor arguments of the stage that has them, so `MatchCombineStage` requires a sampleset and has no cache, and only `ScoreStage` has an imputation floor or a shard.

Genuinely shared behaviour — `load_scores`, `import_match_tables`, `scratch_dir` — sits in a small module of functions, not a base class.

`ScorePipeline` is deleted.

**A stage base class or `Protocol` was rejected.** The four stages do not vary behind a common interface: nothing in this codebase chooses between them at runtime. Nextflow chooses, by invoking a different subcommand. One implementation per slot is a hypothetical seam, not a real one. Worse, a base class would have to hold whatever the four stages have in common, which is how `ScorePipeline` became a union of four stages' needs in the first place — it would re-create the problem one level down.

**An orchestrator that runs all four stages in one process was rejected.** It contradicts ADR 0001: the fan-out belongs to Nextflow because a dask-distributed scheduler needs inter-node network reachability that TREs often do not allow. An in-process orchestrator would also be a second definition of the stage graph, free to diverge from `main.nf`.

## Consequences

- No member of any stage module raises because of which stage constructed it. This is checked by pinning each class's public member set to what its CLI actually calls.
- Each stage's CLI constructs exactly one stage class.
- The four modules share no inheritance, so a change to one stage cannot silently change another. The cost is that a genuinely cross-stage change has to be made four times — accepted, because the stages have diverged more than they have converged, and the compiler cannot tell you which of four unrelated stages a base-class change breaks.
- `ScorePipeline` is removed from `pgscatalog.calc.__all__`. This is a public API break; alpha2 is pre-release, so it is not deprecated or migrated. Alphas are not recorded in `CHANGELOG.md`, which only carries released versions.
- **The stage graph is still defined twice**: `main.nf` in `pgsc_calc`, which is authoritative, and `run_pipeline` in this repo's integration conftest, which is test scaffolding and deliberately not shipped. What binds them is the artefact contract, and those artefacts are not yet named or schema-checked — `matches.parquet` and `match_table.parquet` were structurally interchangeable, so feeding SCORE the ungated one silently scored accessions that failed `min_overlap`. Since closed: SCORE declares only the columns it reads, so the binder refuses the ungated table. The other three artefacts have distinct schemas and no demonstrated silent-failure mode, and remain unmarked.
- Future architecture reviews should not re-suggest extracting a stage base class unless something second actually varies behind the interface. Four concrete modules over two real seams is the decision, not an oversight.
