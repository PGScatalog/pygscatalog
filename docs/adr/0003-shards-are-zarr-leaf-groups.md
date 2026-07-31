# Shards are zarr leaf groups in a nested sampleset/chromosome/file hierarchy

## Status

Accepted. Supersedes [ADR 0002](0002-shard-axes-and-combine.md).

## Context

[ADR 0002](0002-shard-axes-and-combine.md) scoped shards to the variant axis for alpha2 and deferred sample-axis batching (a cohort split into files of e.g. 100,000 samples each) on two blockers: matching a shard independently cannot distinguish a true non-match from an implied homozygous-reference absence, and `fill_missing_dosages` hits the `n_minimum_impute` floor once a shard holds fewer samples than the whole cohort.

Both blockers turned out to be consequences of an assumption ADR 0002 never examined: that every fan-out stage fans out at the same granularity. They are properties of matching *at the leaf*, not of the sample axis.

The zarr hierarchy is already the shard's identity everywhere it matters. `_weight_matrix` builds it with `concat_ws('/', sampleset, filename)`; `targetvariants` keys on `(variant_id, filename, sampleset)`; `get_variants_from_zarr` already recurses the group tree. What ties the codebase to one file per chromosome is not the storage layout but a hardcoded depth of two — `get_variants_from_zarr` yields `group.basename` rather than a relative path, and `unzip_zarr` reconstructs ancestor `zarr.json` paths with `split("/")[0]` and `[0:2]`.

Sample batches are also a weaker concept than ADR 0002 assumed. A batch is one slice of a **unified callset** chopped up by sample for file-size reasons — not an independently genotyped cohort, which would be a separate sampleset. Nothing therefore needs to group by batch, so batch identity need not be modelled at all.

## Decision

A shard is **one zarr leaf group**: the cached genotypes of a single target genome file, identified by its position in the hierarchy rather than by the input file that produced it. The cache is the portable interface between stages — after LOAD, no stage knows whether the input was VCF or bgen.

The hierarchy gains a chromosome level:

```
/{sampleset}/{chromosome}/{filename}/genotypes
                                    /meta
```

Unsplit whole-genome input uses the sentinel chromosome `ALL`, so depth is uniform and fan-out is consistent for every input shape.

**Stages fan out at different levels of the same tree.** MATCH fans out per *chromosome*, holding every leaf's `meta` for that chromosome; it reads 1D variant metadata only, never the 3D genotype arrays, so it does not scale with sample count. SCORE fans out per *leaf*, where the genotype arrays and the memory pressure actually are. This is what a hierarchy buys over a flat list of files, and it dissolves ADR 0002's first blocker: matching against a chromosome-complete variant list means no shard has a partial view.

A chromosome-level MATCH also has the whole chromosome in view, so the invariants that sample batching requires live there. LOAD is fanned out per shard and can observe none of them:

- leaf-path uniqueness — samplesheet validation, before any work
- per-chromosome sample disjointness — MATCH
- ragged batches under `unobserved` semantics — MATCH (see [ADR 0004](0004-absence-semantics-is-declared.md))
- cross-chromosome equality of the sample union — MATCH_COMBINE (relaxed to a subset for chrY by [ADR 0007](0007-ploidy-is-observed-per-call.md): Y is the only chromosome a person can wholly lack, so a male-only chrY leaf is legitimate)

These are scoped to chromosomes **present in the cache**, not chromosomes named in the samplesheet: `cleanup_empty_zarr_group` deletes a group whose index query returned nothing, which happens legitimately whenever the scoring files cover no variants on a chromosome.

Batches are **emergent, not modelled**: a batch is simply one of the files under a chromosome. No samplesheet batch column, no batch identifiers, no declared disjointness. Nothing groups by batch, because `min_overlap` is pooled per `(sampleset, accession)` — a technical batch is a slice of one callset, so a match rate that differs between batches is a data-quality alarm, not grounds for scoring some samples and not others.

For alpha2 the sample axis is supported **only for rectangular batches** (every leaf under a chromosome holding the same variants). Ragged batches require homozygous-reference absence semantics, which is out of scope — see [ADR 0004](0004-absence-semantics-is-declared.md).

Score combination is unchanged from ADR 0002: union each shard's `score_table` and reduce with `GROUP BY (sampleset, accession, sample_id)` and `SUM`. ADR 0002 justified this as "correct because every shard shares the same samples", which is narrower than the reduction actually is. The real invariant is that **every `(sample, variant)` pair appears in at most one leaf**. Variant-axis shards contribute disjoint variants for a shared sample; sample-axis shards contribute disjoint samples. The reduction is agnostic to the shape of the tree, which is why SCORE_COMBINE needs no change to support batching.

(This originally read "exactly one leaf". [ADR 0007](0007-ploidy-is-observed-per-call.md) weakened it: a woman and a chrY variant are a pair that appears in *no* leaf. The reduction was always this tolerant — she simply receives no chrY contribution and a correspondingly smaller `allele_count` — so only the wording was too strong.)

## Consequences

- Four code sites change, three of which delete a hardcoded depth assumption rather than adding one: `TargetGenome._zarr_group_path` gains a level; `get_variants_from_zarr` yields a relative path instead of `group.basename`; `unzip_zarr` loops over ancestors instead of slicing `[0]` and `[0:2]`; `_weight_matrix`'s `concat_ws` takes three parts. `write_dosage_regions` needs nothing — zarr indexes `"a/b/c"` natively.
- `get_sample_ids`' `NotImplementedError` stops being reachable without being removed. It only fires because one process sees every shard at once; a SCORE task is staged with only its own leaf's subtree, so there is nothing left to disagree with. The uniformity assumption stops being load-bearing rather than being fixed.
- `check_chromosomes`' rejection of more than one file without a declared chromosome is **deleted**, not special-cased. Two whole-genome VCFs are no longer ambiguous — they are two sample batches under `ALL`, and the disjointness invariant catches the genuine error.
- The `n_minimum_impute` floor is not a blocker. Rectangular batches impute within a batch of ~100,000 samples, far above the default of 50. The unconditional call in `fill_missing_dosages` (raising even when nothing is missing) remains a pre-existing bug worth a one-line guard.
- **Leaf basename collision is a new silent-corruption risk.** `filename` is `pathlib.Path(target_path).name`, so directories are stripped. Users will organise batched inputs as `chr1/batch1.vcf.gz`, `chr2/batch1.vcf.gz`; nesting keeps those apart, but two batches of the *same* chromosome sharing a basename still collide, and because LOAD is fanned out per shard the collision only appears when the per-shard archives are merged, where one leaf silently overwrites the other. Leaf-path uniqueness validation is required, not optional.
- Zarr caches written by earlier alphas use the depth-2 layout and are not readable under the nested hierarchy. This matters in TREs, where re-running LOAD over a biobank is expensive; alpha2 is pre-release, so the caches are discarded rather than migrated.
- Failing an invariant at MATCH_COMBINE means every LOAD task has already run. A pre-flight step reading VCF headers and bgen `.sample` files could fail faster, but it is a Nextflow-side compute optimisation and can be added later — the authoritative check has to exist at MATCH regardless.
