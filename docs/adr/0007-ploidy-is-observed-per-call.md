# Observe ploidy per call and encode hemizygosity with a sentinel

## Status

Accepted. Relaxes one invariant in [ADR 0003](0003-shards-are-zarr-leaf-groups.md). [ADR 0014](0014-declared-sex-is-checked-never-used.md) amends this decision by adding declared sample sex for QC only. Read *Amendment* below before interpreting declared sex as a reversal.

Library checks support the bgen claims. The Consequences record the checks, their results and the two cases initially left untested.

## Context

At this decision, four restrictions exclude X, Y and MT from scoring. Three are deliberate, and one is accidental:

- `check_chromosomes` rejects target chromosomes outside 1–22 (`cli/load.py`).
- `VALID_CHROMOSOMES` limits which scoring-file positions enter the Cache (`lib/scorefile.py:93`).
- The VCF reader constructs a rectangular array with `np.fromiter(gt_tuples, dtype=(np.uint8, 2))`. It catches `ValueError` and raises `NotImplementedError("Sex chromosome genotypes not supported")` (`lib/cache/_vcf.py:126`). However, numpy 2.3.4 broadcasts haploid `(1,)` to `[1, 1]` instead of raising. This silently doubled called haploid genotypes. Only tuples containing `None` reached the replacement exception, so the LOAD gate actually excluded X.
- `vcf_buffer_variants` sorts regions with `key=lambda t: (int(t[0]), t[1])` (`lib/cache/_vcf.py:65`). The VCF path fails on `int("X")` before handling ploidy. Thus, `test_vcf_buffer_sex_variants` expects `pysam.SamtoolsError` while its docstring claims `NotImplementedError`. It passes for the wrong reason.

Two comments already anticipate ploidy support. `_impute.py:51` contains `n_called_alleles = n_called_genotypes * ZARR_PLOIDY  # replace this with a ploidy mask when X / Y / MT support added`. `_pgs.py:304` contains `* ZARR_PLOIDY  # think about X / Y / MT here too`. Both files are in `lib/score`.

Ploidy belongs to a **call**, not a chromosome or sample. One man can be diploid at a pseudoautosomal X variant and hemizygous a few kilobases away. A model keyed only by `(sample, chromosome)` is incorrect at the PAR boundary.

The formats obtain ploidy differently. The reader interface is expected to support plink and other formats later, so these differences matter.

| format | ploidy source | nullploidy determinable? |
| --- | --- | --- |
| VCF | length of the `GT` field | no — `.` is both a nullploid call and a haploid no-call |
| bgen (layout 2) | explicit per-sample ploidy field, exposed as `BgenVar.ploidy` | yes — `ploidy == 0` |
| bgen (layout 1) | none — the format has no per-sample ploidy | no |
| plink | `.fam` sex column plus chromosome code (23/24/25/26) | yes, by construction |

The previous pipeline used plink and handled Y acceptably. plink's `.bed` and `.pgen` have no ploidy dimension: their genotypes are allele counts. plink derives ploidy from external declarations and marks female chrY calls as missing. It uses declared sex rather than detecting nullploidy. This is a valid alternative design, which this ADR rejects below.

## Decision

**A hemizygous effect allele contributes dosage 1, not 2.** For male non-PAR X, PLINK 2 needs explicit `--xchr-model 1` to match this behaviour. Its default `--xchr-model 2` doubles the coding. Other haploid chromosomes use dosages 0–1. See the [PLINK 2 dosage documentation](https://www.cog-genomics.org/plink/2.0/assoc#xchr-model).

**Observe ploidy per call instead of deriving it from declared sex and a PAR map.** The caller already records this decision. A second derivation could disagree with the file, without a principled way to resolve the conflict. The file records pseudoautosomal calls as diploid, so no separate PAR handling is necessary. This avoids maintaining build-specific PAR coordinates.

**Keep the genotype array shape `(n_variants, n_samples, 2)`.** Unused allele slots receive a new sentinel, distinct from the missing-call sentinel:

```
diploid het          [0,   1  ]
diploid missing      [255, 255]
hemizygous ALT       [1,   254]
hemizygous no-call   [255, 254]
nullploid            [254, 254]

dosage       = sum(gt == effect_allele_idx, axis=2)   # 254 matches no allele index
ploidy       = sum(gt != 254, axis=2)
is_missing   = any(gt == 255, axis=2)
```

Derive the ploidy mask requested by `_impute.py:51` instead of storing it. A stored mask would duplicate existing information in another `(n_variants, n_samples)` array per leaf. That array would be roughly as large as the dosage array. Rename `ZARR_PLOIDY` to `ZARR_MAX_PLOIDY`, because the axis describes the widest call rather than every call's ploidy.

bgen supplies ploidy separately through `BgenVar.ploidy`, but this does not justify a separate stored mask. Its handler encodes ploidy in the genotype array. Reading that encoding gives the same answer for all cases, including missing calls. Separate metadata is a format-specific detail, while the Cache deliberately ignores input format.

**A wholly missing call retains the ploidy declared by `GT`.** pysam returns `(None,)` for `.` and `(None, None)` for `./.`. Thus, VCF can represent the table's hemizygous no-call separately from a diploid missing call.

Imputation uses each missing call's ploidy. Storing a man's missing chrX call as diploid imputes and counts two alleles where he can carry one. This changes the score incorrectly, not merely the array shape.

**Reject `.` when every non-missing call at the variant is diploid.** The error identifies the chromosome and position. Such a variant provides no evidence of hemizygosity to justify `.`. A writer that uses `.` for a diploid no-call would otherwise cause imputation at half the correct dosage.

This check is deliberately limited. If any call is genuinely hemizygous, a length-1 `.` is plausible and the reader honours it. A man's no-call at a mixed-sex chrX variant is conformant output. When every call is `.`, no comparison is possible, so the reader also honours the recorded lengths.

A woman's no-call incorrectly written `.` at a mixed-sex chrX variant passes this check. It then receives half the dosage she can carry. Distinguishing her `.` from a man's would require declared sex, which this ADR rejects for deriving ploidy. The alternative would reject every mixed-sex chrX file with missing calls, including most real biobank chrX data. This residual risk is an accepted limitation.

**Reject ploidy above 2 instead of truncating it.** The genotype axis cannot hold wider calls. bgen permits arbitrary ploidy, and the BGEN specification's reference file contains ploidies 3 and 4. This is therefore a reachable case. Keeping only two alleles would produce a plausible score from an unrepresentable call. Human data should never trigger this limit, which reinforces the need for an error.

**Reject non-autosomal inputs that cannot express ploidy.** Do not assume diploidy. This includes layout-1 bgen, which has no per-sample ploidy field. The same principle applies to mixed-sex chrY in VCF: reject an unexpressible distinction instead of guessing.

**`GenomeFileHandler` returns ploidy encoded in the genotype array.** Each handler privately determines that ploidy:

- `VCFHandler` reads `GT` length.
- `BgenFileHandler` reads the ploidy field.
- A future `PlinkHandler` derives it from `.fam` and chromosome codes.

The Cache retains the format-independent interface from [ADR 0003](0003-shards-are-zarr-leaf-groups.md). This introduces no pipeline-wide sample sex column. The plink handler remains responsible for `.fam`.

`GenomeFileHandler` has two implementations, selected at runtime by `get_file_handler` (`lib/cache/_genomefilehandlers.py:223`). This is an existing interface, not a speculative abstraction. [ADR 0006](0006-one-module-per-stage-no-base-class.md) rejects a stage base class, not this interface.

**Allow chrY to contain a subset of the sample union required by ADR 0003.** Y is the only chromosome that a person can entirely lack. Everyone has at least one X and has MT. This biological exception permits common male-only chrY files without a nullploid encoding. Every sample in such a file is hemizygous.

Readers that can determine nullploidy can still support mixed-sex chrY files. Correct ADR 0003's SCORE_COMBINE rule from "exactly one leaf" to "at most one leaf" for every `(sample, variant)` pair. No leaf contains the pair formed by a woman and a chrY variant. Reduction already tolerates this absence: she receives no Y contribution and a smaller `allele_count`. Only the earlier wording was too restrictive.

**MT is an explicit exception to observed ploidy.** MT is biologically hemizygous, but callers use `1`, `1/1` or `0/1`. Honouring these formats directly could assign dosage 1 or 2 to the same homoplasmic variant. Two callsets for the same people could then disagree.

Handlers collapse homozygous diploid-coded MT calls to hemizygous calls: `0/0`→`0` and `1/1`→`1`. **A heterozygous MT call is a hard error.** It represents heteroplasmy, a mixture within one sample, for which polygenic scoring has no defined dosage. Assigning 1, 0 or missing would invent a convention.

X ploidy differs between people, and the caller records that difference. MT ploidy instead follows mitochondrial biology. This distinction justifies the exception.

## Amendment: check declared sex without using it for scoring

[ADR 0014](0014-declared-sex-is-checked-never-used.md) requires `--sex_map` when a scoring file contains sex chromosome variants. This does not reverse the earlier statement that no pipeline-wide sex column is introduced for deriving ploidy.

Declared sex never determines ploidy or supplies ploidy that the format cannot express. It never affects dosage or scores. The rejection of declared sex plus PAR coordinates remains unchanged, and the repository contains no PAR coordinates. Mixed-sex chrY VCF limitations and a woman's no-call incorrectly written `.` remain accepted. The `_vcf.py` comment and `test_vcf.py` docstring that cite this rejection remain correct.

ADR 0014 adds a check in the opposite direction. A declaration could produce a second answer that disagrees with the file. When the declaration serves only as a check, that disagreement becomes the result. ADR 0014 stops the run instead of selecting either answer. The genotype file remains the sole authority for scoring.

## Considered and rejected

**This ADR rejects declared sex plus a PAR coordinate map**, as in plink. VCF and bgen already record ploidy. Deriving it again discards information and can contradict the input. It would also require a samplesheet column for formats that do not need one and maintenance of genome-build-specific PAR coordinates. However, this remains correct inside a future `PlinkHandler`, because plink provides no alternative.

**This ADR rejects a stored ploidy mask array.** `_impute.py:51` suggests one, but it duplicates the genotype array's information in every leaf.

**This ADR rejects duplicating an allele as `[1, 1]` and halving later.** It still requires a mask identifying samples to halve. Losing that mask would silently double every man's X contribution, adding another failure mode to the stored-mask design.

**This ADR rejects a ragged trailing axis.** It represents variable ploidy directly but breaks uniform slicing through `da.from_zarr(..., chunks=(CHUNK, None, None))[target_idx]` in `_dosage.py:181`.

**This ADR rejects inferring nullploidy from patterns of `.`.** A woman has `.` at every chrY variant, while no-calls can be sporadic. This distinction is incidental. It fails for shards with few Y variants and can silently produce an incorrect inference.

**This ADR rejects deferring MT.** Once heteroplasmy raises an error, MT is the simplest chromosome in this design. It is uniformly hemizygous, without mixed ploidy within variants or nullploidy.

## Consequences

### Dosage and overlap

`sentinel_to_nan` is the highest-risk line in this change (`lib/score/_dosage.py:300`). It maps the missing sentinel to NaN. `calculate_effect_allele_dosage` then marks the whole call as NaN through `da.any(nan_mask, axis=2)` (`:280`). If `254` follows that route, every hemizygous call becomes missing and receives imputation. This could silently replace every man's entire X dosage with plausible numbers. A regression test must detect it.

Enabling X/Y/MT cannot decrease the set of accessions passing `min_overlap`. `load_scoring_files` does not filter chromosomes (`lib/scorefile.py:103`). `score_log_table` uses `score_variant_table LEFT JOIN allele_match_table` (`lib/score/_matchlog.py:81`). Thus, X/Y/MT variants already contribute to the denominator as `unmatched`. Existing passes remain passes, while some failures can become passes.

Passing accessions with X variants can receive different scores because those variants now contribute. This corrects scoring rather than introducing a regression.

`allele_count` needs no schema change for haploidy. It is already a per-sample vector because calls can be missing (`lib/score/_pgs.py:304`). Ploidy adds another reason for counts to differ. `n_matched` remains per-accession and genome-wide (`_pgs.py:238`, retained through `any_value` at `score_combine.py:77`). The number of matched variants remains independent of the sample.

`TargetVariants` validation must accept the new sentinel. Its existing check is `np.isin(self._genotypes, [0, 1, MISSING_GENOTYPE_SENTINEL_VALUE])` (`lib/cache/targetvariants.py:78`).

### bgen evidence

The bgen hard-call selector can no longer use only `variant.is_phased`. A biallelic site with ploidy 1 has two probabilities in both phased and unphased representations. With ploidy 2, those representations have four and three probabilities respectively. Select using `(phased, ploidy)` instead.

Haploid biallelic hard-calling is one `argmax` over two columns in either representation. A file mixing phased and unphased variants established that `is_phased` belongs to each variant. Therefore, evaluate the selector per variant.

Library checks established that `BgenVar.ploidy` is a per-sample `uint8` array containing the true ploidy. Missing samples retain correct ploidy. `variant.probabilities` is rectangular at the widest sample's width. NaN padding occurs only in trailing positions, so it remains distinguishable from real probabilities.

The pinned `bgen==1.9.3` and current `1.9.11` produce byte-identical behaviour. Ploidy support therefore needs no version increase. The library exposes no missingness flag. Detect missing calls when the entire probability row is NaN, which is safe because padding is trailing-only.

Do not hardcode the probability dtype. Unphased probabilities use `float32`. Phased probabilities use `float32` for constant ploidy but `float64` for ragged ploidy.

A ploidy-0 sample and a missing sample have identical all-NaN probability rows. Only `BgenVar.ploidy` distinguishes them. This permits `BgenFileHandler` to support mixed-sex chrY where `VCFHandler` cannot. The handler must read the ploidy accessor rather than infer ploidy from probabilities.

The initial layout-1 rejection evidence came from library source that fills a constant ploidy, not from running a layout-1 file. This was weaker evidence than the other checks. The format contains no ploidy field, so a file-based test would add evidence rather than determine the conclusion. The *As built* section records the later test.

The BGEN specification's reference file contains ploidies 3 and 4. This establishes that the ploidy limit is reachable.

**Correction: phased ploidy-0 bgen is readable and now tested.** The earlier record called this case untestable because "the library's writer refuses to create one". The pinned `bgen==1.9.3` does not refuse it. `BgenWriter.add_variant(..., ploidy=[0, 1, 1, 1], phased=True)` writes successfully.

The reader returns `ploidy == 0` and an all-NaN probability row, as for unphased data. `parse_target_variants` reads the file correctly. It encodes nullploidy before consulting phasing, so a ploidy-0 sample never reaches the `(phased, ploidy)` selector. This explains why phasing cannot affect the result. Implementation of 0025 exposed the correction, which received an independent check.

### Cache and sample boundaries

`add_missing_positions_to_lists` creates padding with `shape=(n_samples, 2)` and the missing sentinel (`lib/cache/targetvariants.py:265`). These rows imply diploidy even on haploid chromosomes. They remain harmless because padded positions have no REF. `update_match_table` filters with `WHERE target.ref IS NOT NULL`.

Thus, padding never receives a `target_row_nr` and never enters dosage or ploidy slices. Missingness provides a second defence: an unmatched position contributes no alleles. Uniformly diploid-shaped padding also cannot trigger the mixed-ploidy guard.

Delete the `NotImplementedError` from `check_chromosomes` (`cli/load.py`) and from `parse_target_variants` (`lib/cache/_vcf.py`). Neither has a line number here because neither remains in the implementation. Rewrite `test_vcf_buffer_sex_variants`, which expected the wrong exception.

This ADR does not support users who divide mixed-sex chrY VCFs into shards. Male-only files work. Mixed-sex support requires a reader that can determine nullploidy, which VCF cannot. As in [ADR 0004](0004-absence-semantics-is-declared.md), the design states its limit instead of fabricating genotypes.

A future plausibility check could detect inconsistent ploidy declarations. This parallels ADR 0004's proposed `homref` check. A variant hemizygous throughout one shard but diploid throughout a sibling shard in the same sampleset indicates inconsistent calls, not biology. This check is not currently worth implementing.

### As built

This section records findings from implementation. They do not change the decisions above.

**The hemizygous no-call encoding `[255, 254]` arrived late because of fixture assumptions.** Initially, `vcf_get_genotypes` stored all wholly missing calls as `[255, 255]`, regardless of `GT` length. A man's missing chrX call therefore became diploid. Under [ADR 0009](0009-pooled-allele-frequency-for-imputation.md), it imputed `2 × freq` for a one-allele call and counted 2 in `allele_count`. This silently inflated the score.

The initial explanation claimed that VCF uses `.` for both haploid and diploid no-calls. That claim is false for pysam, which returns `(None,)` and `(None, None)` respectively. The repository's fixtures used `sample["GT"] = None`, which makes pysam emit `.` even on autosomes. Honouring that length would have made the fixture's missing calls hemizygous, so the fixtures distorted the design.

The three affected fixtures now declare their intended ploidy. bgen was never affected because `BgenVar.ploidy` retains correct ploidy for missing calls.

**The first validation rule rejected a no-call whenever any diploid call shared its variant.** This incorrectly rejected a man's conformant no-call in mixed-sex chrX data. It therefore rejected most real biobank chrX data instead of detecting malformed output.

Restricting the rule to "every non-missing call is diploid" accepts the common valid case. It also accepts a woman's incorrectly written `.`, as documented above. Tests cover both outcomes. A `.` among uniformly diploid calls stops LOAD. A `.` at a mixed-sex chrX variant loads as `[255, 254]`.

**Nullploidy required no reader change.** `bgen_to_hardcall` initialises rows with the no-allele sentinel. It writes `255` only into slots that the sample possesses. Its missing-call branch therefore already produces `[254, 254]` for ploidy 0.

The regression needed to establish the source of that distinction. It checks two identical, wholly NaN probability rows. Swapping only their `BgenVar.ploidy` values must swap the encoded calls.

**Nullploidy exposed a crash when `allele_count == 0`.** A woman scored on a chrY-only accession caused division by zero in `score_avg`. DuckDB produced NaN, and SCORE_COMBINE failed while casting to `DECIMAL(18,6)`. An average over no alleles is undefined.

The generated expression is now `score / NULLIF(allele_count, 0)`, which returns NULL instead of a fabricated zero. This was the first code consequence of [ADR 0003](0003-shards-are-zarr-leaf-groups.md)'s change to "at most one leaf".

**Both readers share one MT helper.** It collapses agreeing diploid-coded pairs to hemizygous calls in place. It rejects disagreeing pairs with an error naming heteroplasmy. A shared helper prevents the formats from implementing different interpretations of MT ploidy.

**The chrY sample-union exception is deliberately asymmetric.** Y no longer adds an expectation row, so its leaf may contain a subset of the sampleset. However, every file still contributes to the union. A Y file cannot introduce a sample absent from every other chromosome. That input still fails, identifying the autosome that lacks the sample. The exception therefore permits male-only files without permitting unrelated sample lists.

**A run established that SCORE_COMBINE needed no change.** The run used an autosomal leaf with 3202 samples and a male-only Y leaf with 1599. It produced exactly 3202 output rows. Women's totals used only the autosomal leaf.

The implementation rejects non-autosomal layout-1 bgen, and a test checks this behaviour. The test uses a manually constructed layout-1 file because the ecosystem has no writer. The reader returns constant ploidy 2, as predicted. This resolves the initial evidence limitation recorded above.
