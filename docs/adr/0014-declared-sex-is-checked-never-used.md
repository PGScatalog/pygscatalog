# Declared sex is checked against observed ploidy, never used to derive it

## Status

Accepted. Amends [ADR 0007](0007-ploidy-is-observed-per-call.md), which rejected a sample sex declaration outright; it stays rejected as an *input to scoring*, and is admitted here as an input to a check.

## Context

[ADR 0007](0007-ploidy-is-observed-per-call.md) observes ploidy per call and states plainly that "no pipeline-wide sample sex column is introduced". That was the right call for scoring and it is untouched here. But it leaves a real gap: nothing in the pipeline ever notices that the sample labelled female is hemizygous across chrX. A sample swap, a mislabelled cohort, or a sex map keyed on the wrong ID space all score cleanly and silently, and X variants now contribute to those scores.

The users who care most about this run in clinical and biobank settings, where an unverified sex label is not an acceptable input to a reported score.

Two facts about this codebase shape every decision below.

**The Cache holds only the positions the scoring files asked for.** That is deliberate (see `CONTEXT.md`, *Cache*) and it is narrower than the source files. So the material available to any sex check is however many X variants a PGS happens to carry — often a handful, frequently none. Loading more X variants to strengthen the check would break what the Cache is.

**Ploidy is already stored per call.** `ADR 0007`'s encoding puts the no-allele sentinel in unused slots, so `(gt != NO_ALLELE).sum(axis=2)` is the ploidy of every call in the shard. Whether a sample is hemizygous on X is a read, not an estimate.

## Decision

**A declared sex is QC input only.** It never derives ploidy, never fills a ploidy a format could not express, and never reaches a dosage or a score. ADR 0007's rejection of "declared sex plus a PAR coordinate map" stands in full, and no PAR coordinates enter the repository.

This is what keeps the two decisions coherent. A declaration that can be *used* is a second answer free to disagree with the file, with no principled way to resolve the conflict — 0007's argument, still sound. A declaration that can only be *checked* has the opposite property: disagreement is the point of it.

**The check is an observed hemizygous fraction, not plink's F coefficient.** For each sample, over the informative X variants in its shard: the fraction of its called genotypes that are hemizygous. Declared male expects 1, declared female expects 0, and the sample is classified on the **majority**.

F was rejected on four counts, the first of which is disqualifying on its own:

- F would be estimated from whatever X variants the scoring files asked for. plink's `--check-sex` defaults assume thousands of LD-pruned markers; a dozen produces a confident-looking number that is noise. The Cache cannot supply more without ceasing to be the Cache.
- F needs allele frequencies, which under [ADR 0009](0009-pooled-allele-frequency-for-imputation.md) are pooled across ploidy — a sex-blended input to a test for sex.
- F needs non-PAR, which needs the build-specific PAR map ADR 0007 declined to maintain.
- F exists because plink's `.bed` has no ploidy dimension at all and heterozygosity is the only signal available. That limitation is not ours.

**PAR is excluded without coordinates, by requiring hemizygosity to be in evidence.** A variant is **informative** when at least one call in the shard is hemizygous. A variant diploid in every sample is excluded, because it is either pseudoautosomal or a cohort with no hemizygous samples in it and nothing distinguishes those two without coordinates.

Excluding them is not a refinement. A scoring file whose X variants are mostly pseudoautosomal makes a man diploid at most of them, landing his fraction near 0 and reading him as female. Without this rule the check is wrong rather than weak.

The consequence is that an all-female shard has no informative variants and cannot be checked. That is correct: it is genuinely indistinguishable from an all-PAR variant set, and reporting either answer would be inventing one.

**Below ten informative variants the check does not conclude, and warns loudly.** Hemizygosity is deterministic rather than estimated, so the floor is not sizing for statistical power — it is tolerance to bad calls. At one informative variant a single miswritten no-call inverts the fraction and stops a run over a sample that is fine. Ten with a majority rule tolerates four.

Not concluding is a warning, not a failure: the alternative is refusing to score every cohort whose PGS lacks X coverage, which is most of them, and would make the mandate self-defeating.

**A mismatch is a hard exception, and `0` is the override.** Where a sample is declared male or female and the majority of its informative calls says otherwise, LOAD raises. Users exempt individual samples by declaring them `0` (unknown).

This reverses an earlier position in this design session and the reason it was reversed should be recorded, because the argument against it is real. Sex chromosome aneuploidy is not rare — Klinefelter (XXY) reads as diploid X, Turner (XO) reads as hemizygous — so a biobank-scale cohort contains correctly declared, correctly genotyped samples that will mismatch. UK Biobank publishes field 31 and field 22001 side by side rather than reconciling them for exactly this reason.

Report-only was rejected anyway, because a warning in a log next to a successful exit code is not a quality check in the setting this feature exists for. Stopping and requiring an explicit per-sample exemption makes the decision a human's, records it in an input file, and keeps it out of the pipeline's hands. The cost is that a legitimate aneuploid sample requires a user to declare it unknown — which is a fair description of what the pipeline knows about it.

**The map is mandatory when a scoring file carries sex chromosome variants**, and is supplied only as `--sex_map`. Keyed on the scoring files rather than on the target genomes because the scoring files are what decide whether the Cache will hold anything to check, and because `Scorefiles.get_unique_positions` (`lib/scorefile.py:70`) can answer it during argument validation, before any genome is read. MT does not trigger it — MT is hemizygous in everyone and says nothing about sex.

**Only the map is read.** No VCF header field records sex. A bgen `.sample` file's `sex` column is deliberately ignored, even where one exists (`tests/data/bgen/phased/tiny1000G.sample` has one). Matching a column by name and reading it under an assumed `0`/`1`/`2` coding is two stacked inferences, and nothing inferred should sit upstream of an exception that stops a run. A mis-read coding would invert every sample and fail the whole cohort, which the user would reasonably read as corrupt genotypes rather than a bad guess by the parser.

This binds the pipeline, not the test suite. The integration fixtures *do* read that column, to write the map the pipeline is then handed. That is the point of them: a map derived from the same genotypes the check reads would make the passing case a tautology, and the `.sample` column is the only declaration in this repository that is independent of the calls.

**LOAD throws after writing the zip, not when the mismatch is found.** The sex map influences nothing that gets cached, so a cache written before the exception is valid and complete. The user fixes the map, re-runs with `--zarr_zip_file`, and `get_unique_positions` returns nothing (`cli/load.py:313`) — the retry costs a check rather than a re-load.

The consequence for the implementation is that the calls are read *before* `shutil.rmtree` deletes the directory store, and the verdict is held until after `zip_zarr`. The check reads the store rather than the `TargetVariants` produced while caching, because on a re-run nothing is cached and that loop never runs.

This looks backwards, which is why it is written down: LOAD deliberately finishes work it is about to reject. Throwing at the point of discovery would make each iteration of an override loop re-read chrX, over a biobank hours at a time, against a workflow this ADR has *designed* to iterate.

**The exception names every mismatching sample at once.** Failing on the first would make a cohort with four hundred mismatches four hundred round trips.

**Map coverage is asymmetric.** A genotyped sample missing from the map is a hard error, because a mandatory declaration that is silently incomplete is not a declaration. A sample in the map that this shard does not hold is fine and expected: sample batches see a subset of the cohort ([ADR 0003](0003-shards-are-zarr-leaf-groups.md)), so one cohort-wide map is the normal way to run this, and rejecting extra rows would make the map a function of how the cohort was chopped up.

## Considered and rejected

**Declared sex as a ploidy source of last resort**, for formats that cannot express it — layout-1 bgen, mixed-sex chrY VCF. Rejected: it spends ADR 0007's strongest guarantee to buy back cases 0007 already documents refusing as a deliberate limit, and it drags the PAR coordinate map back in.

**plink's F coefficient.** Four reasons above; the marker count is the disqualifying one.

**Report-only, never fatal.** The honest option for aneuploidy, and rejected because a clean exit code that quietly means "unverified" is the failure mode this feature exists to remove. The `0` exemption is what makes stopping tolerable.

**A new sex-check stage.** Unnecessary: `CONTEXT.md` guarantees every `(sample, variant)` pair appears in at most one shard, and multiple files under one chromosome are sample batches with disjoint samples, so one chrX shard holds every cached X call for its samples. The check is complete inside a single LOAD and needs no combine. [ADR 0006](0006-one-module-per-stage-no-base-class.md) is not engaged.

**Using chrY as a second signal.** Only bgen can express the nullploidy that would make a female's chrY informative ([ADR 0007](0007-ploidy-is-observed-per-call.md)), so it would make the check's strength format-dependent — a thing to explain forever. X carries it.

**Separating "declared unknown" from "deliberately exempted" into distinct codes.** Rejected for the MVP: `0` is plink and UK Biobank convention and users already know it. The overload is mitigated by reporting, not by a new code.

## Consequences

- **The escape hatch can disable the entire check, and only reporting will show it.** A user facing four hundred hard failures will set the column to `0` and re-run, and the pipeline then exits cleanly having verified almost nothing — worse than report-only ever was, because the exit code now asserts something. Every run must therefore report how many samples were exempted by declaration alongside how many were checked. A run that checked 3 of 500,000 samples has to say so in the same breath as it says it passed.
- **The fraction must be NULL, never `0.0`, when there are no informative variants.** `0/0` computed naively is `0.0`, which reads as diploid, which reads as female-consistent: every woman passes silently and every man is flagged, on a check that never ran. This is the same class of error as `sentinel_to_nan` in ADR 0007 — an arithmetic accident that produces a plausible answer.
- **A scoring file with X variants and a target with no chrX loaded skips the check with a warning.** The map was still mandatory, so the user supplied a file that was never read. Accepted as the price of enforcing the mandate up front rather than deciding requirements from data not yet read. The same applies to a scoring file carrying only chrY variants: the mandate fires, the check has no X to run on.
- **`pgsc_calc run` forwards `--sex_map`** alongside `--bgen_sample_file` (`cli/run.py:96`), because without it `run` cannot score any X-carrying PGS once the mandate lands. Each target is a separate LOAD with its own cache directory, so each writes its own report.
- **Existing chrX integration tests need a sex map**: `test_haploid_chrx.py`, `test_pooled_imputation.py`, and the chrX cases in `test_bgen.py`. `test_male_only_chry.py` needs one too — its scoring file carries chrY, and the mandate keys on X *or* Y — so those tests supply a map and then take the skip-and-warn path. This is the accepted cost of the mandate.
- The flag is `--sex_map`, underscored, matching `--bgen_sample_file`, `--absence_semantics`, `--score_paths` and the rest of `add_load_input_args`.
- A per-sample TSV (`sample_id`, `declared`, `hemizygous_fraction`, `n_informative`, `outcome`) is written whether or not the run throws. It lands at `cache_dir/sex_check.tsv`, beside the archive: LOAD is the one stage with no `--out_dir`, and adding one for a single file would be new surface for Nextflow to thread through. Under Nextflow, LOAD is fanned out and each task reports only its own batch; because batches are disjoint the union is the full picture, but no single log holds it.
- Two outcomes for the MVP — consistent and mismatch. An intermediate fraction (mosaicism, or a broken file) collapses into mismatch rather than getting its own category; the fraction is in the TSV, so it remains diagnosable. Worth splitting later, since it is not actually a finding about sex.
- **A fraction of exactly `0.5` is a mismatch**, not a pass and not an abstention. "Classify on the majority" has no answer there, and it is reachable — ten informative variants, five hemizygous. It is treated as contradictory evidence rather than the absent evidence the floor of ten covers: half the calls say male and half say female, so the sample or the file is wrong either way, and `0` is available to a user who disagrees.
- **A hemizygous no-call (`[255, 254]`) makes a variant informative but appears in neither half of a sample's own fraction.** The file declared ploidy one there, which is the signal informativeness asks about; the fraction's denominator is *called* genotypes, so the sample contributes nothing to its own ratio at that variant.
- **The map needs a header naming `sample_id` and `sex`, and is refused otherwise.** Detecting an optional header would mean deciding from the data whether line 1 is a declaration, which is the same stacked inference this ADR refuses when it declines to read the bgen `sex` column.
- **Coverage is enforced on every target once the mandate fires**, autosomes included, rather than only where the check will run. A map keyed on the wrong sample ID space is one of the three failures this feature exists to catch, and enforcing everywhere makes it fail at the first target instead of whenever chrX happens to be dispatched.
- The whole implementation is one module, `lib/sexcheck.py`, called from `load_cli`. No new stage and no new CLI subcommand — see the rejection above.
- The check reads ploidy the readers already produce, so it adds no format-specific code and inherits ADR 0007's guarantees on both paths. It also inherits 0007's known VCF weakness: a woman's no-call miswritten `.` at a mixed-sex chrX variant is stored hemizygous, so a cohort with enough of them could push a fraction across the majority line. The floor of ten is the only thing standing against that, and it is a real limit rather than a theoretical one.

## Review corrections

The floor and reported `n_informative` apply to each sample's called informative X variants. Missing calls cannot supply the evidence required to make a conclusion. Whole-genome (`ALL`) shards select X rows from cached chromosome metadata; other chromosomes do not contribute to this check.
