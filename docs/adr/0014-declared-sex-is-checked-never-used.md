# Check declared sex against observed ploidy without deriving ploidy from it

## Status

Accepted. This amends [ADR 0007](0007-ploidy-is-observed-per-call.md), which rejected
a sample sex declaration. The declaration remains excluded from scoring inputs.
This decision permits it only as input to a check.

## Context

[ADR 0007](0007-ploidy-is-observed-per-call.md) observes ploidy per call and states
that "no pipeline-wide sample sex column is introduced". This remains the correct
scoring policy. However, nothing detects a sample labelled female that is
hemizygous across chrX.

Several errors can therefore produce scores without a warning:

- A sample swap.
- A mislabelled cohort.
- A sex map that uses the wrong sample identifiers.

X variants now contribute to those scores. Clinical and biobank users need to
check declared sex before accepting a reported score.

Two properties of the codebase constrain the check.

**The Cache holds only positions requested by scoring files.** This is deliberate
(see `CONTEXT.md`, *Cache*). The Cache contains fewer positions than the source
files. A sex check therefore sees only the X variants present in the PGS. There
are often few such variants, and frequently none. Loading additional X variants
would violate the Cache contract.

**The Cache already stores ploidy per call.** ADR 0007 places the no-allele
sentinel in unused slots. Therefore, `(gt != NO_ALLELE).sum(axis=2)` gives each
call's ploidy. The check reads hemizygosity on X directly, rather than estimating
it.

## Decision

**Declared sex is QC input only.** It never determines ploidy or supplies ploidy
that an input format cannot express. It never contributes to dosage or score.
ADR 0007 still rejects declared sex combined with a PAR coordinate map. No PAR
coordinates enter the repository.

Using a declaration to determine ploidy would create a second answer that could
conflict with the file. There is no principled way to resolve that conflict, as
ADR 0007 explains. Checking the declaration has a different purpose. The check is
supposed to detect disagreement.

**The check uses an observed hemizygous fraction, not plink's F coefficient.**
For each sample, it considers informative X variants in that sample's shard. The
fraction counts hemizygous called genotypes divided by all called genotypes at
those variants. Declared male expects 1, and declared female expects 0. The check
classifies the sample by the **majority**.

F was rejected for four reasons. The first is sufficient by itself.

- F would use only X variants requested by the scoring files. Plink's
  `--check-sex` defaults assume thousands of LD-pruned markers. A dozen markers
  produces an apparently precise number that is noise. Supplying more markers
  would violate the Cache contract.
- F needs allele frequencies.
  [ADR 0009](0009-pooled-allele-frequency-for-imputation.md) pools these across
  ploidy. A sex check would therefore use frequencies pooled across sexes.
- F needs non-PAR variants. This requires the build-specific PAR map that ADR
  0007 rejected.
- Plink's `.bed` has no ploidy dimension, so heterozygosity is its only signal.
  This explains its use of F. Our cache does not have that limitation.

**Evidence of hemizygosity excludes PAR without coordinates.** A variant is
**informative** if at least one call in the shard is hemizygous. The check
excludes a variant that is diploid in every sample. Such a variant could be
pseudoautosomal, or the cohort could lack hemizygous samples. Without coordinates,
these cases cannot be distinguished.

This exclusion is necessary for correctness. If most X variants in a scoring
file are pseudoautosomal, a man's calls are mostly diploid. His fraction would
approach 0, and the check would classify him as female. Without the exclusion,
the check would be wrong rather than merely weak.

An all-female shard therefore has no informative variants and cannot support the
check. It is indistinguishable from an all-PAR variant set. Reporting either
classification would invent evidence.

**Below ten informative variants, the check warns and does not conclude.**
Hemizygosity is deterministic, so the floor protects against incorrect calls
rather than setting statistical power. With one informative variant, one
miswritten no-call can invert the fraction and stop a valid sample's run. Ten
variants with a majority rule tolerate four incorrect calls.

Failure to conclude produces a warning, not an exception. Otherwise, most cohorts
would fail because their PGS lacks X coverage. That would defeat the check's
purpose.

**A mismatch raises an exception, and `0` provides an override.** LOAD raises
when a sample's declared male or female sex contradicts its majority of
informative calls. Users exempt individual samples by declaring them `0`
(unknown).

This reverses an earlier design choice despite a valid objection. Sex chromosome
aneuploidy is not rare. Klinefelter (XXY) produces diploid X calls, while Turner
(XO) produces hemizygous calls. Biobank-scale cohorts therefore contain correctly
declared and correctly genotyped samples that will mismatch. For this reason, UK
Biobank publishes field 31 and field 22001 separately rather than reconciling
them.

A report-only check was still rejected. A warning beside a successful exit code
is insufficient quality control for these users. Stopping requires a human to
record each exemption in an input file. The pipeline does not decide the
exemption. A legitimate aneuploid sample must consequently be declared unknown,
which accurately describes the pipeline's knowledge.

**A scoring file with sex chromosome variants requires `--sex_map`.** This is
the only input route for the map. The requirement depends on scoring files,
because they determine which positions the Cache contains.
`Scorefiles.get_unique_positions` (`lib/scorefile.py:70`) can detect the
requirement during argument checks, before reading any genome. MT does not
trigger it. MT is hemizygous in everyone and supplies no information about sex.

**The check reads only the map.** No VCF header field records sex. The pipeline
deliberately ignores a bgen `.sample` file's `sex` column, even when it exists.
`tests/data/bgen/phased/tiny1000G.sample` contains such a column.

Inferring a column's meaning from its name and assuming its `0`/`1`/`2` coding
would combine two assumptions. Neither belongs before an exception that stops a
run. Incorrect coding could invert every sample and fail the whole cohort. A
user could reasonably interpret that as corrupt genotypes rather than an
incorrect parser assumption.

This restriction applies to the pipeline, not its tests. Integration fixtures
read the `.sample` column and create the map supplied to the pipeline. The
column is this repository's only sex declaration independent of genotype calls.
A map derived from the same calls under check would make a passing test a
tautology.

**LOAD raises after writing the zip, not when it first detects a mismatch.**
The sex map changes nothing in the cache. The cache remains valid and complete
when the exception occurs. The user can correct the map and repeat LOAD with
`--zarr_zip_file`. Then `get_unique_positions` returns nothing
(`cli/load.py:313`), so the retry repeats only the check.

The implementation reads calls before `shutil.rmtree` deletes the directory
store. It retains the verdict until after `zip_zarr`. The check reads the store,
not the `TargetVariants` created during caching. A retry caches nothing, so that
loop does not run.

LOAD deliberately completes work before rejecting the declaration. An earlier
exception would force every exemption retry to read chrX again. This could cost
hours per retry in a biobank. The workflow deliberately supports repeated
checks.

**The exception names every mismatching sample together.** Reporting only the
first mismatch would require four hundred retries for four hundred mismatches.

**Map coverage is asymmetric.** A genotyped sample missing from the map causes
an error. A mandatory declaration cannot silently omit samples. Extra map
samples absent from this shard are allowed and expected. Sample batches contain
only part of the cohort ([ADR 0003](0003-shards-are-zarr-leaf-groups.md)). A
cohort-wide map must therefore remain valid regardless of how the cohort is
partitioned.

## Considered and rejected

**Declared sex as a fallback source of ploidy.** Some formats cannot express the
required ploidy, including layout-1 bgen and mixed-sex chrY VCF. ADR 0007 already
accepts those limits. Using declared sex would violate its strongest guarantee
and require the rejected PAR coordinate map.

**Plink's F coefficient.** The four reasons appear above. The insufficient marker
count is enough to reject it.

**Report-only checks without exceptions.** This would accommodate aneuploidy,
but a successful exit could silently mean "unverified". Preventing that outcome
is the feature's purpose. The `0` exemption makes stopping acceptable.

**A separate sex-check stage.** `CONTEXT.md` guarantees that each
`(sample, variant)` pair appears in at most one shard. Files within one
chromosome represent sample batches with disjoint samples. Each chrX shard
therefore contains every cached X call for its samples. One LOAD can complete
the check without a combine step. No new stage is needed under
[ADR 0006](0006-one-module-per-stage-no-base-class.md).

**ChrY as a second signal.** Only bgen can express the nullploidy that would make
a female's chrY informative ([ADR 0007](0007-ploidy-is-observed-per-call.md)).
Using it would make check strength depend on input format and require a lasting
explanation. X supplies the signal instead.

**Separate codes for unknown sex and deliberate exemptions.** The MVP retains
`0`, which users know from plink and UK Biobank. Reporting explains the two
meanings without another code.

## Consequences

- **The override can disable the whole check, and only reporting reveals this.**
  A user facing four hundred failures will set the column to `0` and retry. The
  pipeline then exits successfully while checking almost nothing. This is worse
  than report-only behaviour because the exit code implies a successful check.
  Every run must report both exempted and checked sample counts. A run checking
  3 of 500,000 samples must state that beside its pass result.
- **With no informative variants, the fraction must be NULL, never `0.0`.**
  Naive handling of `0/0` produces `0.0`, which appears diploid and consistent
  with declared female sex. Every woman would pass and every man would fail a
  check that never ran. This resembles the `sentinel_to_nan` error in ADR 0007.
  An arithmetic error creates a plausible answer.
- **X score variants without loaded chrX cause a warning and no check.** The map
  remains mandatory even though the check never reads it. This is the accepted
  cost of enforcing requirements before reading data. The same applies when a
  scoring file contains only chrY variants. The requirement applies, but no X
  variants support the check.
- **`pgsc_calc run` forwards `--sex_map`** with `--bgen_sample_file`
  (`cli/run.py:96`). Otherwise, the new requirement would prevent `run` from
  scoring any PGS with X variants. Each target uses a separate LOAD and cache
  directory. Each LOAD therefore writes its own report.
- **Existing chrX integration tests need a sex map.** These include
  `test_haploid_chrx.py`, `test_pooled_imputation.py` and the chrX cases in
  `test_bgen.py`. `test_male_only_chry.py` also needs a map because its scoring
  file contains chrY. The requirement applies to X or Y. These chrY tests supply
  the map and expect a warning without a check. This is an accepted cost.
- The flag uses underscores: `--sex_map`. This matches `--bgen_sample_file`,
  `--absence_semantics`, `--score_paths` and the other `add_load_input_args` flags.
- LOAD writes a per-sample TSV whether or not it raises an exception. Its columns
  are `sample_id`, `declared`, `hemizygous_fraction`, `n_informative` and `outcome`.
  The file is `cache_dir/sex_check.tsv`, beside the archive. LOAD alone has no
  `--out_dir`. Adding that argument for one file would require another Nextflow
  interface.

  Nextflow runs LOAD separately for each batch. Each task reports only its own
  samples. The batches are disjoint, so their reports together describe the
  cohort. No single log contains the complete result.
- The MVP has two outcomes: consistent and mismatch. An intermediate fraction,
  from mosaicism or a broken file, counts as a mismatch. The TSV retains the
  fraction for diagnosis. A separate category could be useful later because
  this is not actually a finding about sex.
- **A fraction of exactly `0.5` is a mismatch, not a pass or abstention.** Ten
  informative variants with five hemizygous calls can produce this result. A
  majority rule cannot classify it. Half the calls indicate male and half
  indicate female, so the sample or file is wrong. This is contradictory
  evidence, unlike insufficient evidence below the floor of ten. A disagreeing
  user can apply `0`.
- **A hemizygous no-call (`[255, 254]`) makes a variant informative.** Its
  declared ploidy is one, which supplies the required evidence. However, this
  call contributes to neither part of its sample's fraction. The denominator
  counts only called genotypes.
- **The map must have a header naming `sample_id` and `sex`.** The check refuses
  a map without it. Detecting an optional header would require inferring whether
  line 1 declares a sample. This repeats the assumptions rejected for the bgen
  `sex` column.
- **Once required, map coverage applies to every target, including autosomes.**
  It does not depend on whether that target supports the check. Incorrect sample
  identifiers are one of the three errors this feature must detect. Checking
  coverage everywhere exposes that error at the first target, rather than when
  chrX happens to run.
- The implementation occupies one module, `lib/sexcheck.py`, called from
  `load_cli`. It adds no stage or CLI subcommand, as explained above.
- The check reads ploidy produced by existing readers. It adds no format-specific
  code and inherits ADR 0007's guarantees on both input paths. It also inherits
  ADR 0007's known VCF weakness.

  A woman's no-call miswritten as `.` at a mixed-sex chrX variant becomes
  hemizygous. Enough such calls could move the fraction across the majority
  threshold. The floor of ten is the only protection. This is an actual limit,
  not a theoretical concern.
