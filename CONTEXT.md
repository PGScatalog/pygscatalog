# pygscatalog

Python packages (`pgscatalog.core`, `pgscatalog.match`, `pgscatalog.calc`, `pgscatalog.validate`) used by the [PGS Catalog Calculator](https://github.com/PGScatalog/pgsc_calc) (`pgsc_calc`, a Nextflow workflow) to compute polygenic scores. `pgscatalog.calc` is mid-rewrite: a `1.0.0-alpha.1` Python/dask/zarr/duckdb reimplementation (`v3`) is replacing the plink2/Nextflow-based scoring pipeline used by `pgsc_calc` v2 (`release/v2` branch, bugfix-only).

## Language

### Scoring pipeline (`pgscatalog.calc` v3)

**Split-apply-combine (scoring step)**:
The pattern used to scale PGS calculation across a genotype dataset too large to score in one process. SPLIT partitions target genotype data into independent shards (see Shard); the rest of the pipeline runs as four stages: MATCH (per chromosome, over variant metadata only, producing partial match counts and the chromosome's Assessed variants) → MATCH_COMBINE (pools match counts across a sampleset, decides which accessions pass `min_overlap`) → SCORE (per shard, only for accessions that passed) → SCORE_COMBINE (sums partial per-shard scores into one final score per sample per accession). Stages fan out at different levels of the Cache hierarchy on purpose, which is what lets MATCH see a chromosome whole while SCORE stays bounded to one shard's genotypes. Every fan-out stage is implemented as separate Nextflow tasks rather than handled inside a single Python process with a distributed dask cluster — Nextflow's executor plugins (SLURM, cloud batch, etc.) already solve cross-node portability, which matters because every Trusted Research Environment (TRE) configures this differently, and a dask scheduler/worker cluster would need cross-node network reachability that TREs often don't allow. See [ADR 0001](docs/adr/0001-nextflow-fanout-for-split-apply.md), [ADR 0003](docs/adr/0003-shards-are-zarr-leaf-groups.md).
_Avoid_: chromosome split (too narrow — see Shard), two-stage split-apply-combine (it's four stages: MATCH/MATCH_COMBINE gate SCORE/SCORE_COMBINE)

**Shard**:
The cached genotypes of one target genome file — a leaf of the Cache hierarchy, and the unit SPLIT produces and SCORE fans out over. Identified by its place in that hierarchy rather than by the input file that produced it, so a shard is a property of the cache, not of the input. Shards vary along two axes at once: the **variant axis** (which chromosome) and the **sample axis** (which Sample batch). Stages deliberately fan out at different levels of the same hierarchy — MATCH per chromosome, over every batch's variant metadata at once; SCORE per shard, where the genotype arrays are. SCORE_COMBINE sums partial per-shard results grouped by sample ID, valid because every (sample, variant) pair appears in exactly one shard — not because shards share samples. See [ADR 0003](docs/adr/0003-shards-are-zarr-leaf-groups.md).
_Avoid_: chromosome, chunk, file, sample batch (a grouping of samples across shards — see Sample batch)

**Match rate threshold (`min_overlap`)**:
The minimum fraction of a scoring file's variants that must be found in the target genome before an accession's score is calculated (default 75%), used to guard against silently publishing an inaccurate reproduction of a PGS Catalog score. Pooled per (sampleset, accession) — never decided per shard, and never per sample batch: a technical sample batch is a slice of one callset, so a match rate that differs between batches is a data-quality alarm rather than a reason to score some samples and not others. Whether an absent variant counts against the rate depends on Absence semantics.
_Avoid_: overlap threshold, match threshold (use the canonical parameter name `min_overlap`), per-batch match rate (not a gate — see Sample batch)

**Assessed variants**:
The set of variants a callset actually interrogated for a chromosome, whether or not any individual file records them. This is the set MATCH compares scoring files against, and the distinction that Absence semantics turns on: a variant absent from a shard but assessed is homozygous reference, while a variant absent and never assessed is simply unknown. Defined by its role, not its provenance — currently it can only be the union of what the sampleset's shards were observed to contain. Note that an external population reference (e.g. 1000 Genomes, HGDP) cannot supply this set: it records what is polymorphic in humans, not what a given callset interrogated, and treating the two as interchangeable would convert honestly unmatched variants into confident homozygous-reference calls.
_Avoid_: sites, union, variant universe

**Absence semantics**:
What it means for a variant to be absent from a shard: `unobserved` (the variant was never assessed — the historical assumption, rooted in imputed SNP arrays) or `homref` (it was assessed and every sample in that shard is homozygous reference — the sparse-called WGS assumption). Not detectable from a target genome file alone: the two produce identical files, so it is declared per sampleset. It is a single fact with several consequences — it decides whether an absent variant counts against `min_overlap`, whether an absent-but-assessed variant is imputed or filled as a known homozygous-reference dosage, and whether ragged sample batches are an error or expected. Only meaningful for VCF input; bgen is dense, so every bgen shard records every assessed variant.
_Avoid_: missingness (overloaded — a genotype call within a file can also be missing, which is a separate thing), sparse

**Sample batch**:
One slice of a **unified callset** that has been chopped up by sample (e.g. 100,000 samples per file) because a single file covering the whole cohort would be unmanageable. Batches are disjoint in samples and together cover the sampleset. Because every batch comes from the same calling process, batches are expected to agree on which variants exist; where they don't (**ragged** batches — typically because monomorphic sites were dropped per batch), the difference reflects Absence semantics, not a real difference in data quality. A batch is not an independently genotyped cohort — samples genotyped on different arrays or called in separate runs are different samplesets, not batches of one.
_Avoid_: sample set (means Sampleset), chunk, shard (a batch is a grouping of samples; a shard is one file — see Shard), sample axis

**Sampleset**:
A human-readable label for a set of target genomes belonging to the same cohort (e.g. `UKBiobank`), spanning one or more shards. The zarr storage hierarchy is `/{sampleset}/{chromosome}/{shard filename}/genotypes`, chosen so the cache is self-describing: the partitioning is legible to a domain expert inspecting the store, without reference to a samplesheet. Unsplit whole-genome input occupies the sentinel chromosome `ALL`, so the depth is the same whatever the input shape.

**Cache**:
The zarr store written by LOAD and read by every later stage. It is the portable interface between them — once genotypes are cached, no downstream stage knows or cares whether the input was VCF or bgen, and the shard is defined over the cache rather than over the input file that produced it. Consequently checks on the shape of the data (see Sample batch, Absence semantics) are made against the cache, which holds only the positions the scoring files asked for. That is narrower than the source files and deliberately so: it is exactly the subset that can affect a score.
_Avoid_: zarr store (implementation), archive
