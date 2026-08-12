# 0038 — A serial run's peak memory is the max over stages, not their sum

**Depends on:** 0036
**Refs:** ADR 0013, ADR 0003.

## Why

The claim that makes a Serial run usable on real data is that peak memory equals the
largest single stage, not the sum of four — SCORE still holds one shard (ADR 0003) and
MATCH still sees one chromosome, so what a serial run gives up is wall-clock and nothing
else. That claim is what lets `CONTEXT.md` say the limit is shard count rather than
dataset size, and it is currently an argument rather than a measurement.

It holds only if each stage releases what it held. Nextflow gets that right for free by
ending the process. One leaked store, one open duckdb connection or one retained dask
graph turns max into sum, and the symptom would be a serial run dying on data that
`pgsc_calc` handles — the exact failure this entry point exists to avoid.

## Scope

One integration test over `grid_targets`: run every stage individually, record peak RSS
per stage, then run `pgsc_calc run` over the same targets and assert its peak stays within
a small factor of the largest individual stage.

`resource.getrusage(RUSAGE_SELF).ru_maxrss` is enough and needs no dependency. It is
monotonic within a process, which is what makes it the right instrument here: it records
the high-water mark, so a stage that failed to release is visible in the total even
though the sampling is coarse.

Pick the factor from a measurement rather than from taste, and leave the measured numbers
in the test as a comment — the test's job is catching a regression into `sum`, not
pinning an exact number.

## Out of scope

Fixing a leak this uncovers. If the assertion fails on first write, that is its own
ticket with its own diagnosis — this one only builds the instrument.

Per-stage memory limits, and anything that makes `run` enforce `--max_memory_gb` itself.

## Acceptance

- The test passes on the fixture data, with the factor and the observed numbers recorded.
- Artificially retaining a stage's cache across the sequence fails it.
