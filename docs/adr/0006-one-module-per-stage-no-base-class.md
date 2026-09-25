# Use one module per stage without a stage base class or an in-process orchestrator

## Status

Accepted

## Context

One `ScorePipeline` class previously drove all four scoring stages.
[ADR 0001](0001-nextflow-fanout-for-split-apply.md) and [ADR 0003](0003-shards-are-zarr-leaf-groups.md) describe those stages.
The class exposed the combined requirements of all four stages.
Only `load_scores` served more than two stages.

Three interface members depended on the stage, but their type signatures did not show these restrictions:

- `shard` raised an exception for three of the four stages.
- MATCH and SCORE read `sampleset` from the cache.
  MATCH_COMBINE and SCORE_COMBINE had no cache, so their constructors needed `sampleset` explicitly.
- `minimum_samples_for_impute` applied only to SCORE.

Callers needed to know the stage before they could know which members were valid.
The class accumulated checks for these conditions instead of expressing them through construction.

Two shared components already had separate modules.
`GenotypeCache` owned the zarr hierarchy.
`ScoreDatabase` owned connections and each task's resource limits.

## Decision

Four concrete modules represent the stages named in `CONTEXT.md`:

- `MatchStage`
- `MatchCombineStage`
- `ScoreStage`
- `ScoreCombineStage`

Each module accepts only its own requirements and uses the two shared components.
Stage-specific values become constructor arguments only for the stage that needs them.
`MatchCombineStage` therefore requires a sampleset and has no cache.
Only `ScoreStage` has an imputation floor or a shard.

A small function module contains shared behaviour: `load_scores`, `import_match_tables` and `scratch_dir`.
It is not a base class.
The design deletes `ScorePipeline`.

**We rejected a stage base class or `Protocol`.**
The four stages do not provide interchangeable implementations of one interface.
Nothing in this codebase selects between them at runtime.
Nextflow selects a stage by invoking its subcommand.

With one implementation for each role, the proposed abstraction has no actual variation to represent.
A base class would also collect whatever the four stages share.
That approach previously turned `ScorePipeline` into the combined interface of four stages.
Repeating it would restore the same problem at another level.

**We rejected an orchestrator that runs all four stages in one process.**
It conflicts with ADR 0001, which assigns task distribution to Nextflow.
A dask-distributed scheduler needs network access between nodes, which TREs often prohibit.
An in-process orchestrator would also define the stage graph a second time.
That definition could diverge from `main.nf`.

## Consequences

No stage module member raises an exception because of the stage that constructed it.
A test checks each class's public members against the members its CLI calls.
Each stage's CLI constructs exactly one stage class.

The four modules share no inheritance.
Changing one stage therefore cannot silently change another through a shared base class.
A genuinely cross-stage change must occur in four places.
This cost is acceptable because the stages differ more than they agree.
The compiler cannot identify which unrelated stage a base-class change would break.

`pgscatalog.calc.__all__` no longer exports `ScorePipeline`.
This breaks the public API.
Alpha2 is pre-release, so the project provides neither deprecation nor migration.
`CHANGELOG.md` records released versions only, not alphas.

**The stage graph still has two definitions.**
`main.nf` in `pgsc_calc` is authoritative.
This repository's integration conftest contains `run_pipeline` as test scaffolding, which the package deliberately does not ship.
The artefact contract connects these definitions.

At the time of this decision, the artefacts lacked explicit type markers and schema checks.
`matches.parquet` and `match_table.parquet` were structurally interchangeable.
SCORE could therefore accept the ungated table and silently score accessions that failed `min_overlap`.

That defect was subsequently closed.
SCORE now declares only the columns it reads, so the binder rejects the ungated table.
The other three artefacts have distinct schemas and no demonstrated silent-failure mode.
They remain unmarked.

Future architecture reviews should propose a stage base class only when a second implementation actually varies behind the interface.
The decision deliberately uses four concrete modules over two shared components.
