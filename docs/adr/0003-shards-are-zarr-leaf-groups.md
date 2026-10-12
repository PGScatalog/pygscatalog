# Shards are zarr leaf groups in a nested sampleset/chromosome/file hierarchy

## Status

Accepted. Supersedes [ADR 0002](0002-shard-axes-and-combine.md).

## Context

[ADR 0002](0002-shard-axes-and-combine.md) limited alpha2 shards to the variant axis and deferred sample-axis batching.
A sample batch can contain, for example, 100,000 samples from a larger cohort.
That ADR identified two blockers:

- Independent shard matching cannot distinguish a genuine non-match from implied homozygous-reference absence.
- `fill_missing_dosages` can reach `n_minimum_impute` when a shard contains fewer samples than the full cohort.

ADR 0002 assumed that every stage distributed work at the same granularity.
Both blockers followed from matching at the leaf level, not from the sample axis itself.

At the time of this decision, the zarr hierarchy already identified shards throughout the relevant code:

- `_weight_matrix` constructed paths with `concat_ws('/', sampleset, filename)`.
- `targetvariants` used `(variant_id, filename, sampleset)` as its key.
- `get_variants_from_zarr` already traversed the group tree recursively.

A hardcoded depth of two, rather than the storage structure, restricted the code to one file per chromosome.
`get_variants_from_zarr` yielded `group.basename` instead of a relative path.
`unzip_zarr` reconstructed ancestor `zarr.json` paths with `split("/")[0]` and `[0:2]`.

ADR 0002 also assumed that batches needed more structure than they require.
A batch contains a sample subset from one **unified callset**, separated to reduce file size.
An independently genotyped cohort belongs to a separate sampleset.
No operation needs to group by batch, so the model needs no batch identity.

## Decision

A shard is **one zarr leaf group** containing the cached genotypes from one target genome file.
Its position in the hierarchy identifies it, rather than the original input file.
The cache provides a portable interface between stages.
After LOAD, no stage knows whether the original input was VCF or bgen.

The hierarchy gains a chromosome level:

```
/{sampleset}/{chromosome}/{filename}/genotypes
                                    /meta
```

An unsplit whole-genome input uses the sentinel chromosome `ALL`.
Every input shape therefore has the same depth and follows the same task-distribution rules.

**Stages distribute tasks at different levels of the same tree.**
MATCH runs per chromosome and holds every leaf's `meta` for that chromosome.
It reads only 1D variant metadata, never the 3D genotype arrays, so its data requirements do not grow with sample count.
SCORE runs per leaf, which contains the genotype arrays and determines memory demand.

The hierarchy allows matching against a complete chromosome variant list.
Unlike a flat file list with leaf-level matching, no shard has only a partial view.
This removes ADR 0002's first blocker.

Chromosome-level MATCH can check the invariants that sample batching requires.
LOAD runs per shard and cannot check them.
The checks run at these stages:

- Samplesheet validation checks leaf-path uniqueness before any work.
- MATCH checks that samples within a chromosome are disjoint across leaves.
- MATCH checks ragged batches under `unobserved` semantics. See [ADR 0004](0004-absence-semantics-is-declared.md).
- MATCH_COMBINE checks equality of sample unions across chromosomes.
  [ADR 0007](0007-ploidy-is-observed-per-call.md) permits a subset for chrY.
  Y is the only chromosome that a person can wholly lack, so a male-only chrY leaf is valid.

These checks cover chromosomes **present in the cache**, not every chromosome in the samplesheet.
`cleanup_empty_zarr_group` removes a group when its index query returns nothing.
This can legitimately occur when scoring files contain no variants for a chromosome.

The model does not declare batches explicitly.
A batch is simply one file under a chromosome.
There is no samplesheet batch column, batch identifier or declaration of disjointness.

Nothing groups by batch because `min_overlap` pools matches per `(sampleset, accession)`.
A technical batch is a subset of one callset.
Different match rates between batches indicate a data-quality problem.
They do not justify scoring some samples while rejecting others.

Alpha2 supports the sample axis **only for rectangular batches**.
Every leaf under a chromosome must contain the same variants.
Ragged batches need homozygous-reference absence semantics, which remain outside scope.
See [ADR 0004](0004-absence-semantics-is-declared.md).

Score combination remains unchanged from ADR 0002.
It unions shard `score_table` outputs and reduces them with `GROUP BY (sampleset, accession, sample_id)` and `SUM`.
ADR 0002 justified this reduction by saying that every shard shares the same samples.
That justification was narrower than necessary.

The actual invariant is that **every `(sample, variant)` pair appears in at most one leaf**.
Variant-axis shards contribute disjoint variants for shared samples.
Sample-axis shards contribute disjoint samples.
The reduction does not depend on the tree's shape, so SCORE_COMBINE needs no change for batching.

The invariant originally required each pair to appear in exactly one leaf.
[ADR 0007](0007-ploidy-is-observed-per-call.md) relaxed this requirement because a woman and a chrY variant can appear together in no leaf.
The reduction already supported that case.
She receives no chrY contribution and has a correspondingly smaller `allele_count`.
Only the original wording was too restrictive.

## Consequences

Four code locations change, and three changes remove a hardcoded depth assumption:

- `TargetGenome._zarr_group_path` gains a level.
- `get_variants_from_zarr` yields a relative path instead of `group.basename`.
- `unzip_zarr` traverses ancestors instead of selecting `[0]` and `[0:2]`.
- `_weight_matrix` passes three parts to `concat_ws`.

`write_dosage_regions` needs no change because zarr directly supports paths such as `"a/b/c"`.

`get_sample_ids` retains its `NotImplementedError`, but the error becomes unreachable.
Previously, one process saw every shard and could find inconsistent sample lists.
Each SCORE task now receives only its own leaf's subtree.
The uniformity assumption no longer affects the result, although the code still contains it.

`check_chromosomes` no longer rejects multiple files without declared chromosomes.
The rule disappears rather than receiving an exception.
Two whole-genome VCFs represent two sample batches under `ALL`.
The disjointness check detects the actual error.

The `n_minimum_impute` floor does not block this design.
Rectangular batches impute within approximately 100,000 samples, well above the default of 50.
At the time of this decision, `fill_missing_dosages` still raised unconditionally when nothing was missing.
That pre-existing bug warranted a one-line guard.

**Leaf basename collisions can silently corrupt data.**
`filename` uses `pathlib.Path(target_path).name`, which removes directories.
The hierarchy keeps `chr1/batch1.vcf.gz` separate from `chr2/batch1.vcf.gz`.
However, two batches of the same chromosome still collide if their basenames match.

Because LOAD runs per shard, the collision appears only when the shard archives merge.
One leaf then silently overwrites the other.
Leaf-path uniqueness validation is mandatory.

Earlier alpha caches use the depth-2 layout, which the nested hierarchy cannot read.
Repeating LOAD across a biobank is expensive in TREs.
However, alpha2 is pre-release, so the project discards these caches instead of migrating them.

Every LOAD task completes before MATCH_COMBINE can report an invariant failure.
An earlier check of VCF headers and bgen `.sample` files could report failures sooner.
That Nextflow optimisation can follow later.
MATCH must still contain the authoritative check.
