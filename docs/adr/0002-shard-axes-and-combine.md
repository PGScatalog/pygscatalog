# Shards are scoped to the variant axis for alpha2; combine by summing per (sampleset, accession, sample_id)

## Status

Superseded by [ADR 0003](0003-shards-are-zarr-leaf-groups.md) and [ADR 0004](0004-absence-semantics-is-declared.md).

Both blockers recorded below turned out to be consequences of fanning every stage out at the same granularity, rather than properties of the sample axis. Retained for the reasoning, not the decision.

## Context

SPLIT previously assumed a shard was always one chromosome (the variant axis: different variants, same sample cohort across shards). We explored generalizing this to also support the sample axis (sample batches — e.g. 100,000 samples per file, same variants across files, disjoint sample sets), since the underlying zarr storage doesn't structurally forbid it — the hierarchy is already keyed per shard file (`sampleset/filename/genotypes`) — and the only obstacle was code that hard-assumed a uniform sample list across a sampleset's shards (`get_sample_ids()` in `lib/score/_dosage.py` raises `NotImplementedError` if two shards disagree; `check_chromosomes()` in `cli/load.py` rejects more than one unsplit-chromosome file as input).

Investigating what sample-axis support would actually require surfaced problems specific to that axis, not the variant axis:

- **Sparse representation breaks matching correctness.** With sparse-called genotype files, a variant absent from a batch's file usually means homozygous-reference in that batch, not "not genotyped." Matching each shard independently can't tell these apart. Pooling match *counts* sampleset-wide (across batches) fixes the `min_overlap` accept/reject decision, but it doesn't fix scoring itself: a shard's `target_row_nr` values are only ever derived from that shard's own local matching (`_load_variants_from_zarr` reads only the current run's local zarr store), so a variant matched only in another batch is never indexed by this shard at all — no crash, but a **silent** per-sample score bias. Samples in a batch missing a matched-elsewhere variant have that variant's contribution dropped entirely (not 0, not imputed, not counted, not logged), so two samples in the same run can have their score computed from a different effective number of variants with no visible signal that it happened.
- **Small batches break imputation unconditionally.** `fill_missing_dosages` (`_dosage.py`) calls `calculate_mean_dosage` regardless of whether the batch's dosage array has any missing values at all, and `calculate_mean_dosage` raises `ValueError` if `n_samples` is below `n_minimum_impute` (default 50) — before checking whether imputation is even needed. Today this is harmless because `n_samples` is always the whole cohort. Sample-batching makes hitting this floor easy and turns it into a hard crash for small batches even when nothing needed imputing.

Both problems are consequences of a shard only having a partial view of the sampleset when split by sample, and both would need a real fix (matching that distinguishes true non-matches from implied-homref absences; imputation that draws on more than one shard's samples) to be safe to ship.

## Decision

For alpha2, a shard is scoped to the **variant axis only**: one target genotype input file, varying by chromosome (or other variant partitioning) but not by sample subset — i.e. generalizing the existing chromosome-shard model, not extending it. This is what alpha2 actually needs: the motivating problem ([pgsc_calc discussion #459](https://github.com/PGScatalog/pgsc_calc/discussions/459)) was variant-count scaling — combining all chromosomes' genotype data into one unified store for ~200 scores / 12M variants — which variant-axis splitting fully addresses on its own. Sample-count scaling (very large cohorts) is a real, separate problem, deferred to a future release alongside proper homref-aware matching.

Because shards stay same-samples-different-variants, `get_sample_ids()`'s existing assumption holds and doesn't need changing — no `NotImplementedError` to remove, no samplesheet schema change needed yet.

Because no single shard has genome-wide visibility, `min_overlap` accession accept/reject still can't be decided per shard. Matching is split into its own stage ahead of scoring: MATCH runs per shard (fanned out like LOAD/APPLY) and produces a partial match table; a MATCH_COMBINE step pools match counts across all of a sampleset's shards, computes the genome-wide match rate per accession, and decides which accessions pass `min_overlap`. Only accessions that pass are fanned out to the scoring (APPLY) stage — this avoids spending APPLY compute on accessions that are going to be rejected anyway. The pipeline is therefore MATCH → MATCH_COMBINE (decide) → SCORE → SCORE_COMBINE (sum), not a single split→apply→combine pass.

Score combination: each shard's `score_table` output (`sampleset`, `accession`, `sample_id`, `n_matched`, `allele_count`, `score`, `dosage_sum`) is unioned across shards and reduced with duckdb via `GROUP BY (sampleset, accession, sample_id)` with `SUM` over the numeric columns — correct here because every shard shares the same samples, so a given sample's partial contributions from each chromosome shard simply add up.

`min_overlap` itself is a legacy concept from imputed SNP array input, where an unmatched variant genuinely means "not genotyped or imputed" — its purpose is to catch a target dataset too sparse to be a faithful reproduction of the scoring file's published PGS Catalog score (e.g. matching 1% of variants isn't a valid replication). That assumption doesn't hold for sparse-called WGS data, which is part of why sample-axis batching (where the gap actually bites) is deferred rather than shipped with a caveat.

## Consequences

- Sample-axis batching (`sample_batch` samplesheet field, disjoint-sample validation, mixed chrom × batch grids) is explicitly out of scope for alpha2 and needs its own future ADR once homref-aware matching and cross-shard-aware imputation are designed.
- Matching becomes a distinct pipeline stage (MATCH / MATCH_COMBINE) that gates SCORE, rather than something computed inline at the start of a single `pgsc_calc score` invocation as it is today.
- The unconditional `n_minimum_impute` check in `calculate_mean_dosage` (raising even when no imputation is needed) is a pre-existing latent bug independent of this ADR's scope decision — low risk today since `n_samples` is always the full cohort, but worth fixing opportunistically since it will matter as soon as sample-axis batching is picked back up.
