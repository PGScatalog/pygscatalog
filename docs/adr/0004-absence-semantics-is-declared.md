# Absence semantics is declared per sampleset; ragged batches are an error, not a prompt to declare homref

## Status

Accepted

## Context

A variant absent from a target genome file means one of two incompatible things:

- **`unobserved`** — it was never assessed. The historical assumption, rooted in imputed SNP array input, and what `min_overlap` is built on: an unmatched variant genuinely was not genotyped or imputed, so a target too sparse to faithfully reproduce a published PGS Catalog score gets rejected.
- **`homref`** — it was assessed and every sample in that file is homozygous reference. The sparse-called WGS assumption, where emitting all-homozygous-reference rows for 100,000 samples is exactly the file size problem the caller was avoiding.

**The two produce byte-identical files.** Nothing in a VCF distinguishes a chopped-and-AC-filtered callset from an independently joint-called one, so absence semantics cannot be detected — only declared. bgen sidesteps this entirely: it is dense, every sample carries a probability at every variant, so a bgen shard records every variant it assessed and ragged bgen is unrepresentable.

Whether shards actually *disagree* on which variants exist is observable from the cache. Why they disagree is not, and the causes point in opposite directions:

- monomorphic sites dropped when chopping a unified callset (`--min-ac 1`) — absence really is homozygous reference
- a sample subset re-imputed or re-filtered for computational reasons — a variant may have failed QC in the smaller batch, so absence is genuinely unobserved
- independent joint calling per batch — a site is emitted only if some sample in that batch is non-reference, so absence is homozygous reference, but sites with no coverage are indistinguishable

## Decision

Absence semantics is **declared once per sampleset**, as a samplesheet column, defaulting to `unobserved`. Not per shard: mixing semantics within a sampleset would mean the shards are not slices of one callset, which makes them separate samplesets. Not inferred from file type: VCF is both.

The declaration is a single fact with several consequences, which is why it is one flag rather than several knobs. It decides whether an absent variant counts against `min_overlap`, whether an absent-but-assessed variant is imputed or filled as a known homozygous-reference dosage, and whether ragged shards are expected or an error.

For alpha2, **only `unobserved` is implemented**, and ragged shards under `unobserved` semantics are a hard error at MATCH.

**The error must not tell the user to declare `homref`.** Ragged shards are ambiguous in cause: the same shape arises from dropping monomorphic sites (where homozygous reference is true) and from re-imputing a sample subset (where it is false and would silently fabricate genotypes). The user has to know how their callset was generated; the pipeline cannot infer it, and steering them toward a mode that matches the symptom would convert an honest error into a wrong answer. The actionable advice is rectangularity — re-export the batches without dropping monomorphic sites, or split on the variant axis instead.

More generally, `homref` is a claim about a callset's data-generating process, not a remedy for an error message. It should never be documented as the fix for ragged batches.

## Consequences

- Users who chopped a callset *and* filtered monomorphic sites are blocked until `homref` lands. They get a clear error rather than a wrong answer, which is the intended trade.
- `homref` is strictly additive and needs no rework of alpha2. When it lands, a shard scoring a variant it does not contain fills a constant: dosage is `2 * (effect_allele_idx == 0)`, since every sample in that shard is homozygous reference and index `0` is the reference allele. `effect_allele_idx` is already carried per zarr group, so no new match-table columns are needed. The fill rows concatenate onto the real rows before the effect-type adjustment, and because the weight matrix is already pivoted per zarr group, ordering the shard's variants present-first keeps the write contiguous.
- Those fill rows must be materialised, not skipped. When the effect allele is alternate the fill is `0` and contributes nothing to `score` or `dosage_sum`, so it looks droppable — but `allele_count` must still increase by 2 and `n_matched` by 1, or the per-sample average-weighted score is wrong. Skipping them reintroduces exactly the silent per-sample bias ADR 0002 warned about.
- **An external population reference cannot supply absence semantics.** A 1000 Genomes / HGDP callset records what is polymorphic in humans, not what a given callset interrogated. Inferring homozygous reference from it would fill confident genotypes at sites the pipeline never looked at, inflating `n_matched` and `allele_count` and driving the match rate toward 100% by construction — defeating the guard `min_overlap` exists to provide. Restricting to common variants bounds this, since a variant common in 1000G/HGDP is assessed by essentially any modern platform, but it does not close it: common variants dropped by imputation INFO-score QC are routine and would be silently converted from honestly unmatched into confidently called.
- A population reference is still useful for the separate problem of **sampleset absence** — a variant in no shard at all, where the truth is genuinely unknown and `2 × AF` is a better guess than the in-sample mean. That is a different fill from shard absence and needs no coordination with sample batching. It will need its own ADR, because backfilling decouples "matched" from "has a dosage": those are the same set today, and keeping `min_overlap` meaningful means counting backfilled variants as unmatched for the gate while using them for the score.
- A plausibility check is available as a future safeguard against a misdeclared `homref`: if a variant is polymorphic with appreciable frequency in one shard, it is implausible that every sample in a sibling shard of 100,000 is homozygous reference. Not worth building until `homref` exists.
