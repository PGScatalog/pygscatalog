# pygscatalog

The Python packages `pgscatalog.core`, `pgscatalog.match`, `pgscatalog.calc` and `pgscatalog.validate` support the
[PGS Catalog Calculator](https://github.com/PGScatalog/pgsc_calc).
The calculator, `pgsc_calc`, is a Nextflow workflow that calculates polygenic scores.

`pgscatalog.calc` is undergoing a rewrite.
The `1.0.0-alpha.1` implementation, called v3, uses Python, dask, zarr and duckdb.
It replaces the plink2/Nextflow scoring pipeline used by pgsc_calc v2.
The `release/v2` branch accepts bug fixes only.

## Language

### Scoring pipeline (`pgscatalog.calc` v3)

**Split-apply-combine (scoring step)**:
A pattern for calculating polygenic scores when one process cannot score the whole genotype dataset.
SPLIT divides target genotypes into independent shards (see Shard).
Four stages then process the data:

1. MATCH reads variant metadata for each chromosome. It produces partial match counts and the chromosome's Assessed variants.
2. MATCH_COMBINE combines match counts across a sampleset. It decides which accessions pass `min_overlap`.
3. SCORE calculates scores for each shard, using only accessions that passed.
4. SCORE_COMBINE sums partial scores into one final score per sample and accession.

Stages operate at different levels of the Cache hierarchy.
MATCH sees a whole chromosome. SCORE holds only one shard's genotypes.
Nextflow dispatches separate tasks for stages that process multiple parts of the data.
One Python process with a distributed dask cluster does not dispatch those tasks.

Nextflow executor plugins already support different computing systems, including SLURM and cloud batch services.
This matters because each Trusted Research Environment (TRE) configures its systems differently.
A dask cluster would need network access between nodes, which TREs often prohibit.

See [ADR 0001](docs/adr/0001-nextflow-fanout-for-split-apply.md) and
[ADR 0003](docs/adr/0003-shards-are-zarr-leaf-groups.md).

_Avoid_:
- **Chromosome split** is too narrow. See Shard.
- **Two-stage split-apply-combine** is incorrect. MATCH and MATCH_COMBINE control which accessions reach SCORE and SCORE_COMBINE.

**Serial run**:
One process that executes all four scoring stages in sequence over every shard, without an executor to dispatch them.
This entry point serves environments where Nextflow cannot run or where separate tasks add unnecessary complexity.
The stages, artefacts and levels of data partitioning remain unchanged.
A serial run passes the same files between stages as Nextflow tasks use.

Peak memory therefore equals the largest single stage, not the sum of all stages.
Serial execution increases elapsed time.
Its practical limit is the number of shards, not the dataset size.

See [ADR 0013](docs/adr/0013-serial-run-is-a-convenience-entry-point.md).

_Avoid_:
- **Local mode** is ambiguous. Nextflow's default executor is already local.
- **Small-data mode** is incorrect. Shard count limits this mode, not sample or variant count.
- **Orchestrator** obscures the decision in ADR 0006. That ADR rejects an orchestrator as the scaling mechanism. A serial run is a convenience entry point.

**Shard**:
The cached genotypes of one target genome file.
A shard is a leaf of the Cache hierarchy, a unit that SPLIT produces, and a unit that SCORE processes.
Its position in the cache identifies it, rather than the source file alone.
Thus, a shard is a property of the cache.

Shards vary along two axes:

- The **variant axis** identifies the chromosome.
- The **sample axis** identifies the Sample batch.

MATCH reads all batches' variant metadata for one chromosome together.
SCORE reads the genotype arrays for one shard.
SCORE_COMBINE groups partial results by sample ID and sums them.
This sum is valid because each sample/variant pair occurs in at most one shard.
It does not require shards to share samples.
Some pairs occur in no shard, such as a woman and a chrY variant. The sum already permits that absence.

See [ADR 0003](docs/adr/0003-shards-are-zarr-leaf-groups.md).

_Avoid_: **chromosome**, **chunk**, **file** and **sample batch** as synonyms.
A Sample batch groups samples across shards.

**Match rate threshold (`min_overlap`)**:
The minimum fraction of scoring-file variants that must match the target genome before the pipeline calculates an accession's score.
The default is 75%.
This threshold helps prevent publication of an inaccurate reproduction of a PGS Catalog score without a warning.

The pipeline combines counts per sampleset and accession.
It never decides acceptance separately for a shard or sample batch.
A technical sample batch is part of one callset.
A difference in match rates between batches indicates a possible data-quality problem, not permission to score only some samples.
Absence semantics determines whether an absent variant reduces the match rate.

_Avoid_:
- **Overlap threshold** and **match threshold**. Use the canonical parameter name `min_overlap`.
- **Per-batch match rate** as an acceptance gate. See Sample batch.

**Assessed variants**:
The variants that a callset assessed for a chromosome, whether or not an individual file records them.
MATCH compares scoring files against this set.
Absence semantics distinguishes an assessed variant absent from a shard from a variant that the callset never assessed.
The first indicates a homozygous-reference call. The second is unknown.

The set is defined by its role, not its source.
Currently, only the union of variants observed in the sampleset's shards can provide this set.
An external population reference, such as 1000 Genomes or HGDP, cannot provide it.
Such a reference records human variation, not which variants a particular callset assessed.
Substituting that reference would incorrectly convert unmatched variants into homozygous-reference calls.

_Avoid_: **sites**, **union** and **variant universe** as synonyms.

**Ambiguous variant**:
A score variant whose alleles do not resolve its strand when compared with a target variant.
A complemented reading also matches an allele in the target.
A palindromic SNP is the usual example: an `A/T` score variant still reads `A/T` on the opposite strand.

Ambiguity also occurs without an other allele when the target contains both the effect allele and its complement.
Only a Multiallelic variant can provide that combination.
Ambiguity belongs to the pairing of a score variant and target variant, not to either variant alone.
Other score variants at the same palindromic target can remain resolvable.
The pipeline excludes ambiguous variants by default.
Scoring one requires a strand assumption, and an incorrect assumption counts dosage for the wrong allele.

_Avoid_:
- **Palindromic SNP** as a synonym. It excludes other ambiguous pairings.
- **Ambiguous site**. Ambiguity is not a property of the target alone.
- **Strand flip**. That term describes a resolved complemented reading.

**Strand flip**:
A match in which the score and target report the same allele from opposite strands.
The score variant matches its reverse complement instead of its original sequence.
Both strands read 5'→3'. A reverse complement reverses base order and complements each base.
For one base, reversing the order has no effect. For longer sequences, a plain complement is incorrect.

An indel has no strand-flipped reading in this matching model.
Its shared anchor base depends on position and remains leftmost on either strand.
The opposite-strand record therefore has a different position and different alleles.
An operation on the two available allele strings cannot recover that record.

One matching reading resolves the strand.
A second matching reading leaves it unresolved and creates an Ambiguous variant.

See [ADR 0011](docs/adr/0011-strand-flips-are-reverse-complements-and-optional.md).

_Avoid_:
- **Complement** as a synonym. A plain complement can be a different sequence, and matching it is a defect.
- **Flipped variant**. The variant does not change. Its reading changes.

**Match collision**:
Two score variants from one accession that resolve to the same allele of the same target variant.
Two weights then compete for one dosage.
The scoring rows can differ in their alleles, so the Catalog correctly treats them as distinct variants.
The collision becomes visible only after matching.

Best-match selection handles the opposite relationship: several candidate matches for one score variant.
In a collision, each score variant already has one best match.
The pipeline excludes colliding rows instead of choosing a weight.
Choosing a weight would require an assumption about the score author's intention.

_Avoid_:
- **Duplicate variant**. Upstream `is_duplicated` uses the allele pair and correctly does not identify these as duplicates.
- **Not-best match**. That term describes the opposite relationship.
- **Duplicate position**. Two score variants at one position collide only when they resolve to the same allele.

**Multiallelic variant**:
A target variant with more than one alternative allele.
Allele position identifies each allele: reference first, then alternatives in the order recorded by the caller.
Genotype calls refer to these positions instead of allele sequences.
Removing an allele would therefore change the indices of later alleles.

The cache retains alleles that no score can match, including the spanning-deletion `*` and symbolic alleles such as `<NON_REF>`.
Matching skips those alleles without removing them.
A score variant matches at most one alternative allele.
The pipeline counts dosage for the matched allele itself.
A strand flip changes which allele matches, not how dosage is counted.

_Avoid_:
- **Multi-allele site** as a synonym.
- **Exploded variant** for v3. v2's `pgscatalog.match` splits multiallelic variants into biallelic equivalents for plink2. v3 does not.

**Ploidy**:
The number of alleles a sample carries at one variant.
A diploid call has two alleles. A hemizygous call has one.
Examples of hemizygous calls include male non-PAR X, male Y, and MT in all samples.
Ploidy belongs to a call, not to a sample or chromosome.
The same sample can have diploid X calls in PAR and hemizygous X calls a few kilobases away.

The pipeline reads ploidy from each call in the target genome file.
It does not derive ploidy from declared sex and a PAR map.
The caller already recorded that decision.
Deriving ploidy again would create a second result that could disagree with the file.
A hemizygous effect allele contributes dosage 1.
A sample's Allele count therefore depends on its own ploidy, not twice its matched-variant count.

See [ADR 0007](docs/adr/0007-ploidy-is-observed-per-call.md).

_Avoid_:
- **Sex chromosome** as a synonym. X is diploid in half the cohort, and MT is hemizygous without being a sex chromosome.
- **Hemizygous sample**. Hemizygosity belongs to a call.

**Declared sex**:
An external declaration of a sample's sex, distinct from what its genotype calls show.
The pipeline checks this declaration but never uses it to determine scoring ploidy.
Ploidy remains a property observed in each call.
If the declaration and calls disagree, the run stops instead of choosing one source.

The declaration uses the plink and UK Biobank codes: `0` for unknown, `1` for male and `2` for female.
A two-column map from sample ID to sex is the only accepted source.
No VCF header field records sex.
The pipeline deliberately ignores the sex column in a bgen `.sample` file.
Reading a named column with an assumed coding would introduce an inference before a check that can stop the run.

Code `0` also exempts a sample from the Sex check.
It states that sex is unknown and instructs the pipeline to exclude that sample from the check.
Every genotyped sample needs a declaration.
The map may contain more samples than a shard because a sample batch contains only part of the cohort.

_Avoid_:
- **Genetic sex** as an input name. It is the conclusion of a check. UK Biobank fields 31 and 22001 represent this distinction.
- **Gender**. It is a different concept that these files do not record.
- **Sex chromosome** as a synonym. See Ploidy.

**Sex check**:
A comparison of Declared sex with the pattern of a sample's X calls during LOAD.
This check is the only use of declared sex.
It reads Ploidy directly instead of estimating plink's inbreeding coefficient F from heterozygosity.
F compensates for a genotype representation that records allele counts but no ploidy.
The Cache already records ploidy.
Estimating F here would also depend on the small set of X variants requested by the scoring files.

Only informative variants contribute to the check.
An informative variant has evidence of hemizygosity somewhere in the shard.
A variant that is diploid in every sample could be pseudoautosomal or belong to a cohort without hemizygous samples.
The pipeline cannot distinguish those cases without the build-specific coordinates excluded by ADR 0007.
It therefore excludes such variants instead of treating them as evidence of female sex.

The check is required whenever a scoring file contains sex chromosome variants.
Those variants determine whether the Cache can contain data for the check.

_Avoid_:
- **Sex imputation** and **sex inference**. The input declares sex. The check either agrees or stops the run.
- **F coefficient**. The check uses observed ploidy instead.
- **QC** as a synonym. It is broader than this single check and input.

**Allele count**:
The number of alleles a sample carries at the variants matched by an accession.
It is the denominator of that sample's average weighted score.
It is specific to each sample because calls can be missing and Ploidy can differ.
The matched variant count, `n_matched`, is different.
It counts the variants matched by an accession across the genome and is the same for every sample.

_Avoid_: **n_matched** and **allele frequency** as synonyms.

**Absence semantics**:
The declared meaning of a variant's absence from a shard.
The two meanings are:

- `unobserved`: the callset never assessed the variant. This historical assumption comes from imputed SNP arrays.
- `homref`: the callset assessed the variant and every sample in the shard is homozygous reference. This assumption applies to sparse-called WGS.

These cases can produce identical target files, so the pipeline cannot infer the meaning from a file alone.
The user declares it per sampleset.
It determines three behaviours:

1. Whether an absent variant reduces the match rate used by `min_overlap`.
2. Whether an assessed but absent variant receives an imputed dosage or a known homozygous-reference dosage.
3. Whether ragged sample batches indicate an error or an expected representation.

Absence semantics applies only to VCF input.
bgen is dense, so each bgen shard records every assessed variant.
Absence semantics does not describe a variant that a sample cannot carry.
A woman and a chrY variant instead have Ploidy zero.

_Avoid_:
- **Missingness** as a synonym. A genotype call inside a file can also be missing, which is a different condition.
- **Sparse** as a synonym.
- **Nullploid**. It is a property of a call, not of a callset.

**Chromosome**:
A level of the Cache hierarchy and the unit that MATCH processes.
The canonical names are `1`–`22`, `X`, `Y` and `MT`.
The sentinel `ALL` identifies unsplit whole-genome input.
LOAD canonicalises input forms such as `chrX`, `0X`, `23` and `chrM` once.
Later stages therefore do not compare names from different dialects.

The vocabulary excludes `XY`.
Pseudoautosomal variants use `X` because the pipeline observes ploidy rather than deriving it from position.
The separate region name previously indicated where male calls were diploid.

See [ADR 0008](docs/adr/0008-canonical-chromosome-names.md).

_Avoid_:
- **Contig** as a synonym.
- Chromosome names with input-specific prefixes.
- **XY** or **PAR** as a chromosome. PAR is a region, which this model does not name separately.

**Sample batch**:
One part of a unified callset divided by sample, for example into files containing 100,000 samples each.
The division avoids one unmanageable file for the whole cohort.
Batches contain disjoint samples and together cover the sampleset.
Because they come from the same calling process, they should agree on which variants exist.

Ragged batches contain different variants, often because each batch omits monomorphic sites.
This difference reflects Absence semantics rather than a real difference in data quality.
A batch is not an independently genotyped cohort.
Samples genotyped on different arrays or called in separate runs belong to different samplesets, not batches of one sampleset.

_Avoid_:
- **Sample set**. It means Sampleset.
- **Chunk**, **shard** and **sample axis** as synonyms. A batch groups samples across shards. See Shard.

**Sampleset**:
A human-readable label for target genomes from the same cohort, such as `UKBiobank`.
A sampleset can contain one or more shards.
The zarr hierarchy is `/{sampleset}/{chromosome}/{shard filename}/genotypes`.
A domain expert can read the partitioning directly from this hierarchy without a samplesheet.
Unsplit whole-genome input uses chromosome `ALL`, so every input layout has the same hierarchy depth.

**Sampleset allele frequency**:
The proportion of called alleles at a target variant that carry a particular allele, pooled across the sampleset's technical sample batches and observed ploidies.
Repartitioning the same samples into different shards does not change this frequency or the missing dosages imputed from it.

_Avoid_: **Shard frequency** as a synonym; a shard contains only part of the evidence.

**Cache**:
The zarr store that LOAD writes and later stages read.
It provides a portable interface between the stages.
After LOAD creates the cache, later stages do not need to know whether the input was VCF or bgen.
The cache defines the shard, rather than the original input file alone.

Checks of the data layout therefore use the cache. See Sample batch and Absence semantics.
The cache contains only positions requested by scoring files.
This deliberately narrower view contains exactly the subset that can affect a score.

_Avoid_: **zarr store** and **archive** as domain synonyms. They describe implementation or storage rather than the domain concept.
