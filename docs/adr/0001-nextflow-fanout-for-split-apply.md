# Run SPLIT and APPLY as separate Nextflow tasks without a distributed dask cluster

## Status

Accepted

## Context

At the time of this decision, `pgsc_calc` v3's `PGSC_CALC_SCORE` Nextflow process used `label 'process_single'`.
One task collected every shard's cached genotypes and merged them into one zarr store.
It computed all scores as one dask array with the local `scheduler="threads"`.
This scheduler uses several cores on one node.

This design failed at scale.
[pgsc_calc discussion #459](https://github.com/PGScatalog/pgsc_calc/discussions/459) described scoring ~200 scores across ~12M variants.
The unified zarr store reached 811GB and exhausted memory even on an instance with 2TB RAM.

The split-apply-combine scoring design needs to distribute tasks across nodes.
`dask.distributed` can coordinate workers across machines.
However, its scheduler and workers need direct network access to each other.

`pgsc_calc` must run in many Trusted Research Environments (TREs).
Examples include UK Biobank RAP, All of Us, Genomics England and institutional HPC systems.
Their network policies differ and often restrict access.
Nextflow already distributes tasks through executor plugins, including SLURM, LSF, AWS Batch and Google Batch.
This portability is a major reason for using Nextflow in `pgsc_calc`.

## Decision

Nextflow runs SPLIT and APPLY as separate task instances, with one task per shard.
The TRE's configured Nextflow executor dispatches these tasks.
This follows the existing `PGSC_CALC_LOAD` design, which also runs one task per shard.

Dask provides parallel computation within a node through the existing local threaded scheduler.
The design introduces no dask-distributed scheduler or worker cluster.

## Consequences

- Additional Nextflow jobs or nodes provide horizontal scaling in every TRE that Nextflow already supports.
  The design adds no network requirements between nodes.
- COMBINE becomes a separate step.
  Each shard's `PGSC_CALC_SCORE` task no longer produces the final combined output directly.
- Each shard remains limited by the resources of one Nextflow task.
  This decision does not address very large individual shards.
