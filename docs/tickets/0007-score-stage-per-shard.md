# 0007 — SCORE stage, fanned out per shard

**Depends on:** 0001, 0004
**Blocks:** 0008
**Refs:** ADR 0001, ADR 0003

## Why

`PGSC_CALC_SCORE` is `label 'process_single'` today: one task merges every shard's genotypes into a single zarr store and scores it as one dask array. This is the OOM in [pgsc_calc#459](https://github.com/PGScatalog/pgsc_calc/discussions/459) — ~200 scores / ~12M variants produced an 811GB store that failed on a 2TB-RAM instance. Fanning out per shard is the fix, and it is the point of alpha2.

## Scope

- `pgscatalog-calc score` operates on **one shard**: its own leaf's subtree only, staged via the existing partial extraction (`unzip_zarr(group_path=...)`).
- Takes the passing accession list from MATCH_COMBINE; rejected accessions cost no compute.
- Emits a partial `score_table` (`sampleset`, `accession`, `sample_id`, `n_matched`, `allele_count`, `score`, `dosage_sum`) instead of writing final output.
- Dask stays on the local threaded scheduler — intra-node only, no dask-distributed cluster. Horizontal scaling is Nextflow's job, because TREs vary too much in what cross-node networking they permit (ADR 0001).

## Notes

`get_sample_ids`' `NotImplementedError` (`_dosage.py:495`) is not reachable here: it fires only when one process sees every shard at once, and a SCORE task sees one leaf. Do not "fix" it — the assumption simply stops being load-bearing.

`get_n_samples` becomes shard-local, which is correct: the dosage array is `(n_variants_for_this_shard, n_samples_in_this_shard)`.

## Acceptance

- Scoring a sampleset shard-by-shard and combining (0008) reproduces today's single-task scores bit-for-bit on the existing test data.
- Peak memory per task scales with one shard, not the sampleset.
- An accession rejected by `min_overlap` is never read.
