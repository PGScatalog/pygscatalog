# Declare absence semantics per sampleset and reject ragged batches without suggesting homref

## Status

Accepted

## Context

An absent variant in a target genome file has two incompatible interpretations:

- **`unobserved`** means that the callset never assessed the variant.
  This historical assumption comes from imputed SNP array input and underlies `min_overlap`.
  An unmatched variant was genuinely neither genotyped nor imputed.
  The check rejects a target that is too sparse to reproduce a published PGS Catalog score faithfully.
- **`homref`** means that the callset assessed the variant and every sample in the file is homozygous reference.
  Sparse-called WGS uses this assumption to avoid writing all-homozygous-reference rows for 100,000 samples.

**Both interpretations can produce byte-identical files.**
A VCF cannot distinguish a subdivided, AC-filtered callset from an independently joint-called one.
The pipeline therefore cannot detect absence semantics.
The user must declare them.

Dense bgen avoids this ambiguity.
Every sample has a probability at every variant, so a bgen shard records every variant that it assessed.
It cannot represent ragged bgen batches.

The cache shows whether shards disagree about which variants exist.
It does not show why they disagree.
Possible causes imply different absence semantics:

- Removing monomorphic sites when subdividing a unified callset, with `--min-ac 1`, makes absence genuinely homozygous reference.
- Re-imputing or re-filtering a sample subset for computational reasons can remove a variant that failed QC in the smaller batch.
  Its absence is genuinely unobserved.
- Independent joint calling emits a site only when a sample in that batch is non-reference.
  Absence therefore means homozygous reference, but sites without coverage remain indistinguishable.

## Decision

A samplesheet column declares absence semantics **once per sampleset** and defaults to `unobserved`.
Shards cannot declare different semantics within one sampleset.
Such shards would not represent subsets of one callset and would belong to separate samplesets.
File type cannot determine the declaration because VCF supports both interpretations.

One declaration controls several related behaviours:

- Whether an absent variant counts against `min_overlap`.
- Whether the pipeline imputes an absent-but-assessed variant or assigns a known homozygous-reference dosage.
- Whether the declaration permits ragged shards or treats them as an error.

These consequences follow from one fact, so the interface uses one flag rather than separate controls.
Alpha2 implements **only `unobserved`**.
MATCH reports a hard error for ragged shards under that declaration.

**The error must not suggest declaring `homref`.**
Ragged shards can result from removing monomorphic sites, where homozygous-reference absence is correct.
They can also result from re-imputing a subset, where that assumption would silently invent genotypes.
The user must know the process that generated the callset.
The pipeline cannot infer that process from the error.

Suggesting a declaration that fits the symptom could turn a valid error into an incorrect score.
The actionable advice is to make the batches rectangular.
Users can re-export batches without removing monomorphic sites, or split along the variant axis instead.

`homref` describes the process that generated a callset.
It is not a remedy for an error message.
Documentation must never present it as the correction for ragged batches.

## Consequences

Users who subdivided a callset and removed monomorphic sites remain blocked until `homref` support exists.
They receive a clear error instead of an incorrect score.
That trade-off is intentional.

Future `homref` support needs no redesign of alpha2.
For a variant absent from a shard, every sample in that shard is homozygous reference.
The dosage is therefore the constant `2 * (effect_allele_idx == 0)`.
Index `0` denotes the reference allele.

Each zarr group already carries `effect_allele_idx`, so the match table needs no new columns.
The implementation adds fill rows after real rows and before the effect-type adjustment.
The weight matrix already pivots per zarr group.
Placing each shard's present variants first therefore keeps writes contiguous.

The implementation must materialise these fill rows, even when the effect allele is alternate.
An alternate effect gives dosage `0` and adds nothing to `score` or `dosage_sum`.
However, it must still increase `allele_count` by 2 and `n_matched` by 1.
Otherwise, the sample's average-weighted score is incorrect.
Omitting these rows would restore the silent per-sample bias described in ADR 0002.

**An external population reference cannot determine absence semantics.**
A 1000 Genomes or HGDP callset records which variants are polymorphic in humans.
It does not identify the sites assessed by a particular callset.
Inferring homozygous reference would assign confident genotypes to sites that the pipeline never assessed.

That inference would inflate `n_matched` and `allele_count`.
It would also drive the match rate towards 100% by construction and defeat the purpose of `min_overlap`.
Restricting the inference to common variants limits the problem but does not eliminate it.
Essentially every modern platform assesses variants common in 1000G or HGDP.
However, imputation INFO-score QC routinely removes common variants.
The inference would silently convert those genuine non-matches into confident calls.

A population reference remains useful for **sampleset absence**, where no shard contains a variant.
The true genotype is genuinely unknown, and `2 × AF` is a better estimate than the in-sample mean.
This operation differs from filling shard absence and needs no coordination with sample batching.

Sampleset filling needs a separate ADR because it separates “matched” variants from variants with dosage values.
Those sets are currently identical.
To preserve `min_overlap`, the overlap check must treat population-filled variants as unmatched while scoring can still use their dosages.

A future plausibility check could detect some incorrect `homref` declarations.
A variant can be polymorphic at appreciable frequency in one shard.
It is then implausible that all 100,000 samples in a sibling shard are homozygous reference.
This check is not worth implementing before `homref` exists.
