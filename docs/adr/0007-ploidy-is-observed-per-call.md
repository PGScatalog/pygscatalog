# Ploidy is observed per call; hemizygosity is a sentinel in the existing ploidy axis

## Status

Accepted. Relaxes one invariant in [ADR 0003](0003-shards-are-zarr-leaf-groups.md). Amended by [ADR 0014](0014-declared-sex-is-checked-never-used.md), which introduces a sample sex declaration for QC only — see *Amendment* below before reading the rejection of declared sex as a reversal.

The bgen claims below were checked against the library rather than assumed — see the Consequences for what was run, what it showed, and the two cases left untested.

## Context

X, Y and MT are excluded from scoring today in four places, three deliberate and one accidental:

- `check_chromosomes` rejects any target genome outside 1–22 (`cli/load.py:131`)
- `VALID_CHROMOSOMES` restricts which scoring file positions get cached (`lib/scorefile.py:93`)
- the VCF reader builds a rectangular array with `np.fromiter(gt_tuples, dtype=(np.uint8, 2))`, whose `ValueError` is caught and re-raised as `NotImplementedError("Sex chromosome genotypes not supported")` (`lib/cache/_vcf.py:126`). That guard is weaker than it looks: numpy broadcasts a haploid `(1,)` to `[1, 1]` rather than raising (checked on numpy 2.3.4), so a haploid **called** genotype was silently doubled, and only a `None` in the tuple ever reached the re-raise. The LOAD gate, not this, is what kept X out
- `vcf_buffer_variants` sorts its region list with `key=lambda t: (int(t[0]), t[1])` (`lib/cache/_vcf.py:65`), so the VCF path actually dies on `int("X")` *before* reaching the ploidy problem. This is why `test_vcf_buffer_sex_variants` asserts `pysam.SamtoolsError` while its docstring claims `NotImplementedError` — the test passes for the wrong reason.

Two sites already anticipate the work: `n_called_alleles = n_called_genotypes * ZARR_PLOIDY  # replace this with a ploidy mask when X / Y / MT support added` (`lib/score/_impute.py:51`), and `* ZARR_PLOIDY  # think about X / Y / MT here too` (`lib/score/_pgs.py:304`).

Ploidy is not a property of a chromosome or of a sample. It is a property of a **call**: the same man is diploid at a pseudoautosomal X variant and hemizygous a few kilobases outside it. Any model keyed on (sample, chromosome) is wrong at the PAR boundary.

Three formats determine ploidy three different ways, which matters because the reader interface is expected to grow a plink implementation and others after it:

| format | ploidy source | nullploidy determinable? |
| --- | --- | --- |
| VCF | length of the `GT` field | no — `.` is both a nullploid call and a haploid no-call |
| bgen (layout 2) | explicit per-sample ploidy field, exposed as `BgenVar.ploidy` | yes — `ploidy == 0` |
| bgen (layout 1) | none — the format has no per-sample ploidy | no |
| plink | `.fam` sex column plus chromosome code (23/24/25/26) | yes, by construction |

plink is worth stating precisely because the previous pipeline used it and handled Y acceptably. plink's `.bed`/`.pgen` has **no ploidy dimension at all** — genotypes are allele counts. It derives ploidy from declarations outside the genotype data and sets female chrY calls to missing. It did not detect nullploidy; it was told the sex. That is a legitimate design, and it is the alternative rejected below.

## Decision

**A hemizygous effect allele contributes a dosage of 1**, not 2. This matches plink2 `--score`'s default; the doubled coding is `--xchr-model 2`, an opt-in.

**Ploidy is observed per call**, never derived from a declared sample sex plus a pseudoautosomal region map. The caller already made this decision and wrote it down; re-deriving it produces a second answer that is free to disagree with the file, with no principled way to resolve the conflict. It also makes PAR free — a pseudoautosomal variant is diploid in the file because it is diploid — which avoids shipping build-specific PAR coordinates and keeping them correct.

**The genotype array keeps its shape.** It stays `(n_variants, n_samples, 2)`; unused allele slots hold a new sentinel distinct from the missing-call sentinel:

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

The ploidy mask that `_impute.py:51` calls for is **derived, not stored**. A stored mask would be a second `(n_variants, n_samples)` array in every leaf — roughly the size of the dosage array — holding information the genotype array already encodes. `ZARR_PLOIDY` becomes `ZARR_MAX_PLOIDY`, since the axis is now the width of the widest call rather than the ploidy of every call.

Deriving the mask survives contact with bgen, which hands ploidy over on a *separate* channel from the probabilities (`BgenVar.ploidy`) and so might look like an argument for storing it. It is not: the handler encodes that ploidy into the array, and reading it back gives the same answer for every case, including missing calls. The separate channel is a fact about one format, and the Cache is deliberately format-blind.

**A wholly missing call keeps the ploidy its `GT` field declares.** The two spellings are distinguishable — pysam returns `(None,)` for `.` and `(None, None)` for `./.` — so the table's `hemizygous no-call` is reachable from a VCF rather than collapsing into `diploid missing`. It has to be: a missing call is imputed at its own ploidy, so a man's missing chrX call stored as diploid is imputed and counted at two alleles he cannot carry, which is a wrong score rather than a wrong shape.

**Where hemizygosity is nowhere in evidence at a variant, a `.` is rejected rather than read.** Reading ploidy from the `GT` length trusts every upstream writer to be conformant, and a writer that spells a diploid no-call `.` would have those samples imputed at half the dosage they can carry — the same class of error in the other direction, and just as quiet. The trust is therefore checked against the rest of the variant: if a `.` appears and **every** non-missing call at that variant is diploid, the reader raises, naming the chromosome and position. That is the signature of a non-conformant writer, because nothing at that position is hemizygous for the `.` to be.

The check is deliberately no stronger than that. If some call at the variant is genuinely hemizygous, hemizygosity is real there and a length-1 `.` is a plausible reading of it, so the length is taken at face value — and a mixed-sex chrX variant with a man's no-call is ordinary conformant output, not a broken writer. A variant where every call is `.` has nothing to compare against at all and is likewise read as written.

The residual risk is stated rather than removed: at a mixed-sex chrX variant, a woman's no-call miswritten `.` passes this check and is imputed at half the dosage she can carry. Telling her `.` from a man's needs a declared sample sex, which this ADR refuses on principle, so the alternative to accepting it is refusing every mixed-sex chrX file that has any missing call at all — which is most real biobank chrX data. This is an accepted limitation of observing ploidy per call, not an oversight.

**Ploidy above 2 is rejected, never truncated.** The genotype axis is two wide, so a higher-ploidy call cannot be held. bgen permits arbitrary ploidy and the BGEN spec's own reference file exercises 3 and 4, so this is reachable rather than theoretical. Keeping the first two alleles would turn a call the Cache cannot represent into a plausible score. Human data should never trip this, which is exactly why it must be an error rather than a comment.

**Inputs that cannot express ploidy are rejected for non-autosomes rather than assumed diploid.** This covers layout-1 bgen, which has no per-sample ploidy field at all. The rule is the same one applied to mixed-sex chrY in VCF: where a format cannot say, the pipeline refuses instead of guessing.

**`GenomeFileHandler` returns ploidy already encoded in the array.** How a handler knew is private to it: `VCFHandler` reads the `GT` length, `BgenFileHandler` reads the ploidy field, a future `PlinkHandler` derives it from `.fam` and the chromosome code. The Cache therefore stays the format-blind interface [ADR 0003](0003-shards-are-zarr-leaf-groups.md) describes, and no pipeline-wide sample sex column is introduced — plink's `.fam` stays plink's business. `GenomeFileHandler` is a genuine seam: two implementations exist and `get_file_handler` chooses between them at runtime (`lib/cache/_genomefilehandlers.py:223`). [ADR 0006](0006-one-module-per-stage-no-base-class.md)'s rejection of a stage base class does not apply here, and should not be read as forbidding this one.

**chrY relaxes ADR 0003's cross-chromosome sample-union invariant to a subset.** Y is the only chromosome a person can lack entirely — everyone has at least one X, and everyone has MT — so the carve-out is biology rather than convenience. This admits male-only chrY files, which are common and wholly unambiguous: every sample in one is hemizygous, so no nullploid encoding is needed at all. Mixed-sex chrY files remain supportable through the nullploid encoding wherever a reader can determine it. ADR 0003's justification for SCORE_COMBINE is corrected from "every `(sample, variant)` pair appears in **exactly** one leaf" to "**at most** one leaf": a woman and a chrY variant are a pair no leaf holds. The reduction was always this tolerant — she receives no chrY contribution and a correspondingly smaller `allele_count` — so only the wording was too strong.

**MT is a deliberate exception to observed ploidy.** MT is biologically hemizygous, but callers write it as `1`, `1/1`, or `0/1`. Honouring the file would let a formatting convention change a score: the same homoplasmic variant reads as dosage 1 from one caller and 2 from another, so two callsets of the same people disagree. Handlers therefore collapse homozygous diploid-coded MT calls to hemizygous (`0/0`→`0`, `1/1`→`1`, unambiguous), and **a heterozygous MT call is a hard error**. A het MT call is heteroplasmy — a within-sample mixture — and polygenic scoring has no defined dosage semantics for it. Calling it 1, or 0, or missing would invent one.

Unlike X, where ploidy genuinely differs between people and the caller is the authority on which, MT ploidy is a fact about mitochondria. That is what makes this an exception rather than an inconsistency.

## Amendment: a sex declaration exists, and is never used

[ADR 0014](0014-declared-sex-is-checked-never-used.md) adds a mandatory `--sex_map` when a scoring file carries sex chromosome variants. A reader who finds that flag and this ADR's "no pipeline-wide sample sex column is introduced" together should not conclude one of them lost.

Nothing above is reversed. Declared sex never derives ploidy, never fills a ploidy a format could not express, and never reaches a dosage or a score. The rejection of *declared sex plus a PAR coordinate map* below is unchanged, no PAR coordinates are in the repository, and the limitations this ADR accepts on principle — the mixed-sex chrY VCF, the woman's no-call miswritten `.` — are still accepted rather than resolved by the new declaration. The `_vcf.py` comment and `test_vcf.py` docstring that cite this ADR as refusing a declared sample sex remain correct as written.

What 0014 adds is the opposite operation. This ADR's argument against a declaration is that it produces a second answer free to disagree with the file, with no principled way to resolve the conflict. A declaration that can only be *checked* inverts that objection: the disagreement is the output, and 0014 stops the run rather than picking a winner. The genotype file remains the sole authority on how a sample is scored.

## Considered and rejected

**Declared sex plus a PAR coordinate map** (plink's model). Rejected because both VCF and bgen already record ploidy, so this discards real information in favour of a derivation that can contradict it, requires a mandatory samplesheet column for formats that need no such thing, and puts genome-build-specific PAR coordinates into the repository to be maintained. It remains the right mechanism *inside* a `PlinkHandler`, where the format leaves no alternative.

**A stored ploidy mask array.** What `_impute.py:51` suggests, and rejected for size: it duplicates what the genotype array already says, once per leaf.

**Duplicating the allele (`[1, 1]`) and halving downstream.** Requires a mask to know which samples to halve, so it is the stored-mask option with an extra failure mode — losing the mask silently doubles every man's X contribution.

**A ragged trailing axis.** Honest about the data, but breaks the uniform `da.from_zarr(..., chunks=(CHUNK, None, None))[target_idx]` slicing in `_dosage.py:181`.

**Inferring nullploidy from the pattern of `.` calls** — a woman is `.` at *every* chrY variant, a no-call is sporadic. Rejected: it infers a claim from an incidental difference, it breaks when a shard holds few chrY variants, and it fails silently when wrong.

**Deferring MT.** Rejected because MT is the *easiest* chromosome under this design — uniformly hemizygous, no mixed ploidy within a variant, no nullploidy — once the heteroplasmy question is answered with an error rather than a guess.

## Consequences

- **The single most dangerous line in the change is `sentinel_to_nan` (`lib/score/_dosage.py:300`).** It maps the missing sentinel to NaN, and `calculate_effect_allele_dosage` then NaNs the whole call via `da.any(nan_mask, axis=2)` (`:280`). If `254` is allowed down that path, every hemizygous call becomes a missing call and is imputed — every man's entire X chromosome, silently, producing a plausible number. This needs a test that fails loudly, not a comment.
- Enabling X/Y/MT is **monotone on `min_overlap`**. `load_scoring_files` applies no chromosome filter (`lib/scorefile.py:103`) and `score_log_table` is `score_variant_table LEFT JOIN allele_match_table` (`lib/score/_matchlog.py:81`), so X/Y/MT scoring variants are already in the match-rate denominator today, counted `unmatched`. No accession that passes now can fail afterwards; only failures can become passes. Score *values* will change for passing accessions that carry X variants, because those variants begin contributing. That is the correctness improvement, not a regression.
- `allele_count` absorbs haploidy without a schema change. It is already a per-sample vector (`lib/score/_pgs.py:304`) because calls go missing; ploidy is simply a second reason it varies. `n_matched` stays per-accession and genome-wide (`_pgs.py:238`, carried through by `any_value` at `score_combine.py:77`), because "how many variants matched" is still sample-independent.
- `TargetVariants`' genotype validation (`np.isin(self._genotypes, [0, 1, MISSING_GENOTYPE_SENTINEL_VALUE])`, `lib/cache/targetvariants.py:78`) gains the new sentinel.
- The bgen reader can no longer select its hard-call function from `variant.is_phased` alone. For a biallelic site, ploidy 1 yields **two** probabilities under both the phased and the unphased representation, where ploidy 2 yields four and three respectively. The selector becomes `(phased, ploidy)`. Convenient consequence: haploid biallelic hard-calling is a single `argmax` over two columns and is identical for phased and unphased input. `is_phased` is per *variant*, not per file — verified on a file mixing both — so the selector has to be evaluated per variant.
- **Verified against the library rather than assumed.** `BgenVar.ploidy` is a per-sample `uint8` array carrying the true ploidy, and it stays correct for samples whose call is missing. `variant.probabilities` is rectangular at the width of the widest sample, padded with `NaN` in **trailing** positions only, so padding is always distinguishable from a real probability. Behaviour is byte-identical between the pinned `bgen==1.9.3` and the current `1.9.11`, so **no version bump is needed for ploidy support**. Two details worth carrying into the implementation: there is no exposed missingness flag, so a missing call has to be detected as *"the entire probability row is NaN"* (safe, because padding is trailing-only); and the `probabilities` dtype is not stable — unphased is `float32`, phased is `float32` at constant ploidy but `float64` when ploidy is ragged, so it must not be hardcoded.
- **A ploidy-0 sample and a missing sample have identical all-NaN probability rows in bgen**, separable only through `BgenVar.ploidy`. This is the whole reason `BgenFileHandler` can support mixed-sex chrY where `VCFHandler` cannot, and it means the handler must read ploidy from the accessor and never infer it from the probability array.
- Layout-1 bgen falls under the rejection rule above. The evidence is weaker than for the rest of this section: it comes from reading the library's source, which fills a constant, rather than from running a layout-1 file. The conclusion follows from the format having no such field, so a test file would confirm rather than inform it.
- Ploidy 3 and 4 were observed in the BGEN spec's reference file, which is what moves the ploidy bound from theoretical to reachable.
- **Correction — a phased ploidy-0 bgen variant is readable, and is now tested.** This bullet originally recorded the case as untestable because "the library's writer refuses to create one". On the pinned `bgen==1.9.3` it does not refuse: `BgenWriter.add_variant(..., ploidy=[0, 1, 1, 1], phased=True)` writes, and the reader returns `ploidy == 0` with an all-NaN probability row, exactly as in the unphased case. `parse_target_variants` reads such a file correctly. It works because a ploidy-0 sample never reaches the `(phased, ploidy)` selector — it is encoded nullploid before phasing is consulted, which is why phasing cannot matter here. Found while implementing 0025 and confirmed independently.
- **Verified:** `add_missing_positions_to_lists` fabricates padding rows at `shape=(n_samples, 2)` filled with the missing sentinel (`lib/cache/targetvariants.py:265`). On a haploid chromosome those rows imply diploidy for every sample, and they are harmless for a stronger reason than being all-missing: a padded position has no REF, and `update_match_table` filters candidates on `WHERE target.ref IS NOT NULL`, so it never gets a `target_row_nr` and is never sliced into the dosage or ploidy arrays. Being all-missing is the second line of defence — an unmatched position contributes no alleles either way. They also can't trip the mixed-ploidy guard, being uniformly diploid-shaped.
- `check_chromosomes`' `NotImplementedError` (`cli/load.py:131`) and the `NotImplementedError` in `parse_target_variants` (`lib/cache/_vcf.py:126`) are both deleted. `test_vcf_buffer_sex_variants` is rewritten — it asserted the wrong exception.
- Users who chopped a mixed-sex chrY VCF are not served by this ADR. A male-only chrY file works; a mixed-sex one depends on the reader determining nullploidy, which VCF cannot. This is the same trade [ADR 0004](0004-absence-semantics-is-declared.md) makes: a clear limit rather than a fabricated genotype.
- A misdeclared ploidy has a plausibility check available later, mirroring the one ADR 0004 describes for `homref`: a variant hemizygous in every sample of one shard and diploid in every sample of a sibling shard of the same sampleset is a callset inconsistency, not biology. Not worth building yet.

### As built

Recorded after implementation. The decisions above stand; these are what shipping them showed.

- **The encoding table's `hemizygous no-call [255, 254]` shipped late, and the fixtures were why.** `vcf_get_genotypes` first stored *every* wholly missing call as `[255, 255]`, whatever its `GT` length, so a man's missing chrX call was a missing **diploid** call: under [ADR 0009](0009-pooled-allele-frequency-for-imputation.md) it imputed `2 × freq` where he can carry one allele, and counted 2 toward his `allele_count` — an inflated score, silently. The reason first recorded for it — that VCF writes `.` for both a haploid and a diploid no-call — **is false** of pysam, which returns `(None,)` and `(None, None)` respectively. The claim came from this repo's own fixtures, built with `sample["GT"] = None`, which makes pysam emit `.` even on an autosome; honouring length would have made every one of those missing calls look hemizygous, so the fixtures were quietly deciding the design. The three that did this now write the ploidy they mean. bgen was never affected: `BgenVar.ploidy` stays correct for missing calls.
- **The check was first written to fire whenever any diploid call shared the variant, and that was wrong.** It refused a mixed-sex chrX file the moment a man had a no-call — which is conformant output from an ordinary caller, so the rule rejected most real biobank chrX data rather than the broken writer it was aimed at. Narrowing it to "every non-missing call is diploid" costs the woman-miswritten-as-`.` case above and buys back the common one. Both sides are tests: a `.` among uniformly diploid calls stops LOAD, and a `.` at a mixed-sex chrX variant loads and is stored `[255, 254]`.
- **Nullploidy needed no reader change at all.** `bgen_to_hardcall` initialises every row to the no-allele sentinel and writes `255` only into the slots a sample actually has, so ploidy 0 falls out as `[254, 254]` from the missing-call branch. What needed proving was not that the encoding works but that the *distinction* is sourced correctly — the regression test asserts the two probability rows are identical and wholly NaN, and that swapping only the `BgenVar.ploidy` entries swaps the output calls.
- **Nullploidy made `allele_count == 0` reachable, and it crashed.** A woman scored on a chrY-only accession divides by zero in `score_avg`; DuckDB produced NaN and SCORE_COMBINE's export died casting it to `DECIMAL(18,6)`. An average over no alleles is undefined, so the generated column is now `score / NULLIF(allele_count, 0)` — NULL rather than a fabricated zero. This is the first place the "at most one leaf" weakening of [ADR 0003](0003-shards-are-zarr-leaf-groups.md) had a consequence in code rather than in wording.
- **The MT exception is one helper shared by both readers**, not two implementations that can drift. It collapses an agreeing diploid-coded pair to hemizygous in place and raises on a disagreeing one, naming heteroplasmy. Keeping it out of the format-specific readers is what makes "MT ploidy is a fact about mitochondria" a single statement in the code.
- **The chrY sample-union relaxation is asymmetric, deliberately.** Y no longer contributes an expectation row, so a chrY leaf may cover a *subset* of the sampleset — but the union is still built from every file, so a Y file can never introduce a sample no other chromosome has. That still fails, and blames the autosome that lacks the sample. The carve-out admits male-only files without also admitting an unrelated sample list.
- **SCORE_COMBINE needed no change, confirmed by running it** rather than by reading the SQL: a sampleset with an autosomal leaf over 3202 samples and a male-only chrY leaf over 1599 produces exactly 3202 output rows, women summed from the autosomal leaf alone.
- Layout-1 bgen rejection for non-autosomes is implemented and tested against a hand-rolled layout-1 file, which the ecosystem has no writer for. The evidence caveat above is now discharged: such a file reads back with a constant ploidy of 2, as predicted.
