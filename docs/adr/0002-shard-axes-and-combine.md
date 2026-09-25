# Limit alpha2 shards to the variant axis and sum scores per sample and accession

## Status

Superseded by [ADR 0003](0003-shards-are-zarr-leaf-groups.md) and [ADR 0004](0004-absence-semantics-is-declared.md).

Both blockers below resulted from running every stage at the same granularity.
They were not inherent limits of the sample axis.
This document retains the original reasoning, not the superseded decision.
The implementation findings below describe the code examined for that decision.

## Context

SPLIT previously assumed that each shard contained one chromosome.
These variant-axis shards contained different variants for the same sample cohort.
We considered extending shards to the sample axis.
For example, files could each contain 100,000 samples, with the same variants and disjoint sample sets.

The zarr structure already supported this arrangement through the per-file hierarchy `sampleset/filename/genotypes`.
The obstacle was code that required the same sample list across a sampleset's shards.
`get_sample_ids()` in `lib/score/_dosage.py` raised `NotImplementedError` when two shards disagreed.
`check_chromosomes()` in `cli/load.py` rejected more than one file without a chromosome split.

Investigation found two problems specific to sample-axis support.

**Sparse representation could make matching incorrect.**
In sparse-called genotype files, an absent variant usually means that the batch is homozygous reference at that variant.
It does not usually mean that the variant was not genotyped.
Independent shard matching cannot distinguish these cases.

Pooling match counts across batches corrects the sampleset-wide `min_overlap` decision.
It does not correct scoring.
Each shard derives its `target_row_nr` values only from its own matches.
`_load_variants_from_zarr` reads only the current run's local zarr store.
A variant matched only in another batch therefore receives no index in this shard.

This omission causes silent per-sample score bias rather than a crash.
The shard drops the variant's contribution entirely.
It neither assigns zero nor imputes, counts or logs the missing contribution.
Two samples in one run can therefore use different effective variant counts without a visible warning.

**Small batches could fail imputation even without missing calls.**
`fill_missing_dosages` in `_dosage.py` called `calculate_mean_dosage` whether or not the dosage array contained missing values.
`calculate_mean_dosage` raised `ValueError` when `n_samples` was below `n_minimum_impute`, whose default was 50.
It checked this limit before checking whether imputation was necessary.
At the time, `n_samples` always described the whole cohort, so this behaviour caused no problem.
Sample batching could easily reach this limit and crash even when every call was present.

Both problems arose because a sample shard saw only part of the sampleset.
Safe support would need matching that distinguishes genuine non-matches from implied homozygous-reference absence.
It would also need imputation that uses samples from more than one shard.

## Decision

For alpha2, a shard contains one target genotype file and varies only along the **variant axis**.
It can represent a chromosome or another variant partition, but not a sample subset.
This extends the chromosome-shard model without adding another axis.

The motivating problem was variant-count scaling, as described in [pgsc_calc discussion #459](https://github.com/PGScatalog/pgsc_calc/discussions/459).
Combining every chromosome produced one genotype store for ~200 scores and 12M variants.
Variant-axis splitting fully addresses that problem.
Sample-count scaling for very large cohorts is a separate problem.
A future release must address it with matching that understands homozygous-reference absence.

All shards retain the same samples and different variants.
The existing `get_sample_ids()` assumption therefore remains valid.
This decision requires neither removal of `NotImplementedError` nor a samplesheet schema change.

No shard sees the whole genome, so it cannot decide whether an accession passes `min_overlap`.
Matching becomes a separate stage before scoring:

1. MATCH runs per shard, like LOAD and APPLY, and produces a partial match table.
2. MATCH_COMBINE pools match counts across every shard in the sampleset.
3. MATCH_COMBINE calculates each accession's genome-wide match rate and applies `min_overlap`.
4. SCORE, also called APPLY, processes only accepted accessions.
5. SCORE_COMBINE sums the partial scores.

The graph is MATCH → MATCH_COMBINE (decide) → SCORE → SCORE_COMBINE (sum).
It is not one split→apply→combine pass.
Rejecting accessions before APPLY avoids scoring work that the overlap check would discard.

Each shard emits a `score_table` with `sampleset`, `accession`, `sample_id`, `n_matched`, `allele_count`, `score` and `dosage_sum`.
DuckDB unions these tables and groups them with `GROUP BY (sampleset, accession, sample_id)`.
It applies `SUM` to the numeric columns.
This reduction is correct here because every shard contains the same samples.
Each sample's partial chromosome contributions add to its final score.

`min_overlap` originated with imputed SNP array input.
For that input, an unmatched variant genuinely means that the variant was neither genotyped nor imputed.
The check rejects target data too sparse to reproduce the scoring file's published PGS Catalog score.
For example, matching 1% of variants does not constitute a valid replication.

Sparse-called WGS data does not satisfy this absence assumption.
The problem becomes material with sample-axis batching.
This decision therefore defers batching instead of releasing it with a caveat.

## Consequences

- Alpha2 explicitly excludes sample-axis batching.
  This includes a `sample_batch` samplesheet field, checks for disjoint samples and mixed chromosome × batch grids.
  A future ADR must follow designs for homref-aware matching and imputation across shards.
- MATCH and MATCH_COMBINE become separate pipeline stages that gate SCORE.
  Previously, one `pgsc_calc score` invocation computed matches before scoring.
- The unconditional `n_minimum_impute` check was a pre-existing bug, independent of this scope decision.
  It raised even when no imputation was necessary.
  The whole-cohort sample count made its risk low at the time.
  Correcting it opportunistically was still worthwhile because sample-axis batching would expose it.
