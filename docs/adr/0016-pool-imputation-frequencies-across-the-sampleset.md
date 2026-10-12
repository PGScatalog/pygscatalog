# Pool imputation frequencies across the sampleset

## Status

Accepted and implemented. Extends [ADR 0009](0009-pooled-allele-frequency-for-imputation.md): pooling across observed ploidies remains, but the estimation population becomes the sampleset rather than the individual shard.

LOAD records per-allele counts sufficient to derive called-allele totals for each cached target variant in each shard. Before SCORE, these counts are summed across the sampleset for the same variant and allele. The pooled frequency is the summed allele count divided by the summed called-allele total, not an unweighted average of shard frequencies. SCORE uses this shared frequency multiplied by each missing call's observed ploidy.

Sample batches are technical partitions of one callset. Repartitioning the same samples must not change their imputed dosages. Keeping shard-local frequencies would preserve existing results but make those results depend on arbitrary file boundaries and smaller estimation samples. Sample-axis runs may therefore change scores; that is an intended correction.

Population-aware allele frequencies are deferred to a future release. Sampleset pooling can overestimate or underestimate missing dosages for groups whose allele frequencies differ from the pooled estimate; this decision does not claim to resolve population structure.

A call with any missing allele contributes nothing to either the per-allele counts or the called-allele total. For example, `0/.` contributes neither its known reference allele nor a denominator of one. This preserves SCORE's existing whole-call missingness policy. Absent ploidy slots are not missing alleles: a fully called hemizygous genotype contributes one allele, and a nullploid call contributes none.

The imputation sample floor remains configurable with a default of 50, measured as the number of distinct samples in the sampleset, not the number in a shard or the number of fully called samples at each variant. Do not introduce a per-variant called-sample threshold. Variants with no called alleles still have no estimable frequency; the sampleset-size check does not make a zero denominator valid. A run with no missing dosages to impute does not require the sample floor.

If a matched variant requires imputation but has no called alleles anywhere in the sampleset, SCORE must stop with a hard error identifying the variant. Do not drop its contribution or substitute an invented frequency. A shard with entirely missing calls can still use evidence from other shards. Nullploid calls are not missing genotypes and must not acquire an imputed allele contribution.

Compatibility with older caches and cache migration are out of scope for this pre-release change. Do not add a migration path or a fallback that recalculates counts from older caches.

LOAD stores one variable-length, slash-separated count vector per variant in `meta/allele_counts`, ordered REF then ALT as recorded in the cache. The called-allele total is the sum of that vector. This uses the existing one-dimensional metadata interface and supports different numbers of alleles across LOAD batches without padding every variant to the maximum allele count.

MATCH already holds every sample shard for its chromosome. It sums the cached counts per target variant and allele before joining frequencies to score matches, so multiple scores cannot duplicate the evidence. MATCH_COMBINE carries these frequencies through the gated match table and adds the distinct sampleset sample count from its existing sample manifests. SCORE consumes these values. No new pipeline stage or intermediate file is needed.

This decision does not extend the supported absence semantics.
