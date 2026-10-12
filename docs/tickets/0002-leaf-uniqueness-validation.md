# 0002 — Validate leaf-path uniqueness; drop the unsplit-chromosome restriction

**Refs:** ADR 0003

## Why

`filename` is `pathlib.Path(target_path).name` (`targetgenome.py:78`), so directories are stripped. Users organise batched inputs as `chr1/batch1.vcf.gz`, `chr2/batch1.vcf.gz`. Nesting (0001) keeps those apart, but two files of the *same* chromosome sharing a basename still collide — and because LOAD is fanned out per shard, the collision only appears when the per-shard archives are merged, where one leaf silently overwrites the other. No error, no NaN, just a sampleset quietly missing genotypes.

Deriving the group name from the parent directory does not fix this: Nextflow stages `path` inputs flat into the task work directory, so the parent is gone by the time Python sees it (except for s3 URIs passed as strings, which would then behave differently from staged files); `/runA/x/f.vcf` and `/runB/x/f.vcf` still collide; and `/` is a zarr path separator, so a name containing it would create an unintended nesting level.

## Scope

- Validate that `{chromosome}/{filename}` is unique across a sampleset's inputs, before any work runs. Pure path arithmetic — no file access.
- Error names both offending inputs and says to rename or use `stageAs`.
- Delete the `check_chromosomes` rule rejecting more than one file without a declared chromosome (`load.py:99`). Two whole-genome VCFs are no longer ambiguous — they are two sample batches under `ALL`, and per-chromosome sample disjointness (0005) catches the genuine error.

## Acceptance

- Two inputs resolving to the same leaf path fail before LOAD, naming both.
- Same basename under different chromosomes is accepted.
- Two whole-genome VCFs with disjoint samples are accepted and land under `ALL`.
