# 0037 — Delete `run_pipeline`; the integration tests drive `pgsc_calc run`

**Depends on:** 0036
**Refs:** ADR 0013, ADR 0006, ADR 0005.

## Why

ADR 0006 tolerated `run_pipeline` in `tests/integration/conftest.py` as a second
definition of the stage graph *because it is not shipped*. Ticket 0036 ships one. Two
definitions is the cost ADR 0006 warned about, and the fix is that the shipped one is the
one the tests exercise — which is also the dogfooding argument for building it at all.

## Scope

Delete `run_pipeline` from `tests/integration/conftest.py`. Its fifteen call sites invoke
`pgsc_calc run` instead, through the existing `run_cli` helper.

Keep `load_targets`, `load_target` and `run_cli`: about ten tests drive a *single* stage
deliberately (`test_four_stage_pipeline.py` lines 391-506, `test_male_only_chry.py`,
`test_nested_hierarchy.py`), and those are testing a stage, not a pipeline. They stay
as they are.

Tests that assert against the directory tree keep asserting against it — the prefix moves
to `out_dir/stages/`, nothing else changes:

- `test_dosage_array_is_sized_to_one_shard` reads `score/N/scratch`
- the ADR 0005 layout test pins each stage directory's contents
- `test_match_logs_are_written` reads `match_combine/logs`

They must not end up pinning a scratch path, which ADR 0005 declares unstable — the tree
under `stages/` is `run`'s output, and each stage directory below it still holds exactly
its declared artefacts plus `scratch`.

## Out of scope

Rewriting what the fifteen tests assert. This is a prefix change and a call swap; any
test whose assertions need rethinking is a signal that 0036's layout is wrong, not
licence to change the test.

## Acceptance

- `run_pipeline` no longer exists, and no test constructs a stage sequence of its own.
- The full integration suite passes unchanged in what it asserts, including the golden
  scores and the ADR 0005 layout test.
- `grep -rn "scratch" tests/integration` shows the same reaches past a stage interface as
  before, and no new ones.
