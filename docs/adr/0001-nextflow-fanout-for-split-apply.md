# Fan out SPLIT and APPLY as Nextflow tasks, not a distributed dask cluster

## Status

Accepted

## Context

`pgsc_calc` v3's `PGSC_CALC_SCORE` Nextflow process is currently `label 'process_single'`: one task collects every shard's cached genotype data, merges it into a single zarr store, and computes scores over the whole thing as one dask array (local `scheduler="threads"`, i.e. parallel across cores on one node only). This fails at scale — see [pgsc_calc discussion #459](https://github.com/PGScatalog/pgsc_calc/discussions/459), where scoring ~200 scores / ~12M variants produced an 811GB unified zarr store that ran out of memory even on a 2TB-RAM instance.

Reimplementing split-apply-combine for the scoring step needs a way to fan work out horizontally. Dask has a native answer to this (`dask.distributed`, a scheduler coordinating workers across multiple machines), but that requires the scheduler and workers to be directly network-reachable to each other. `pgsc_calc` is designed to run inside many different Trusted Research Environments (TREs) — UK Biobank RAP, All of Us, Genomics England, institutional HPC, etc. — each with different, often locked-down network policies. Nextflow already solves "run N tasks across a cluster" portably via its executor plugins (SLURM, LSF, AWS Batch, Google Batch, ...), which is a large part of why `pgsc_calc` is a Nextflow pipeline in the first place.

## Decision

SPLIT and APPLY are implemented as separate Nextflow task instances, one per shard, dispatched through whichever Nextflow executor the TRE is already configured to use — mirroring how `PGSC_CALC_LOAD` already fans out one task per shard today. Dask stays scoped to intra-node parallelism only (the existing local threaded scheduler); no dask-distributed scheduler/worker cluster is introduced.

## Consequences

- Horizontal scaling means adding more Nextflow-managed jobs/nodes, which works in any TRE Nextflow already supports, without new inter-node network requirements.
- COMBINE becomes its own step, since `PGSC_CALC_SCORE` per shard no longer produces final combined output directly.
- Per-shard computation is still bounded by the resources of a single Nextflow task (vertical scaling within a shard) — very large individual shards are not addressed by this decision.
