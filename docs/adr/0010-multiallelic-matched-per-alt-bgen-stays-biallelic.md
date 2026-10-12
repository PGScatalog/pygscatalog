# Match multiallelic variants per alternative allele. Keep bgen biallelic

## Status

Accepted. Defines how score variants match targets with more than one alternative allele. [ADR 0011](0011-strand-flips-are-reverse-complements-and-optional.md) defines a complemented reading and when the matcher attempts one.

## Context

The recorded v2 design in `pgscatalog.match` supports only the first alternative allele because it must provide plink2 with a biallelic pvar. `annotate_multiallelic` splits `ALT` at commas and explodes the frame (`_match/preprocess.py:42-59`). Matching then seeks the biallelic equivalent of each score variant. This representation loses information and raises a concern about strand flips changing dosage after that split.

v3 avoids this problem. `match_variant` returns an `effect_allele_idx`. `calculate_effect_allele_dosage` counts genotype calls equal to that index (`lib/score/_dosage.py`). The matcher compares allele sequences once. Scoring subsequently uses index arithmetic.

Complementing `A` to `T` does not change an allele's position in the target record. `REFALT_FLIP` and `ALTREF_FLIP` therefore return the same index as their corresponding unflipped matches. A flip can change which allele matches, but cannot change how scoring counts that allele. v3 needs no allele-explosion step and must not add one.

Multiallelic support instead requires changes at the Cache boundary. The bcftools filter must retain these records. bgen needs an explicit policy, and the matcher must handle alleles that are not nucleotide sequences.

The existing `--exclude-types other` filter does not remove the problematic symbolic alleles. A direct bcftools check showed that `G,<NON_REF>` passes unchanged. Thus, `*` or `<NON_REF>` can reach matching in multiallelic and other records.

## Decision

**Read multiallelic records and consider every ALT as a candidate.** Match the score variant against each ALT and select the highest priority. No tie-break policy is necessary because at most one distinct ALT can match a given score variant.

With an other allele, a forward match on `alt_i` and complemented match on `alt_j` require `alt_i = reverse_complement(ref) = alt_j`. Without an other allele, the effect allele equals exactly one distinct ALT.

**Retain unmatchable alleles in their original positions and skip them during matching.** Genotypes refer to allele indices. Removing `*` from `[C,*,G]` changes `G` from index 3 to index 2 and makes its dosages refer to the wrong allele.

The spanning-deletion `*` and symbolic alleles such as `<NON_REF>` remain in the allele list. The matcher skips them without an ambiguity check and returns `NO_MATCH` for them. A score of `A/G` against `REF=A ALTS=[*,G]` matches `G` at index 2 and scores normally. Any other non-nucleotide allele string is corrupt and raises an error.

**Ambiguity depends on the score-target pairing.** A score variant is ambiguous when its complemented reading also matches an allele in the target. This generalises the palindrome check, so no separate palindrome test is necessary. The matcher detects ambiguity from a second matching reading.

Every example below uses single bases. Plain complement and reverse complement therefore produce the same result.

| target | score variant | forward | complemented | verdict |
|---|---|---|---|---|
| `A / T` | `A/T` | idx 0 | idx 1 | ✅ ambiguous — palindrome |
| `A / C` | `A/C` | idx 0 | — | ✅ resolvable |
| `A / C` | `C` (no other allele) | idx 1 | — | ✅ resolvable |
| `A / C,G` | `C` (no other allele) | idx 1 | idx 2 | ✅ ambiguous — only reachable multiallelically |
| `A / T,C` | `C/A` | idx 2 | — | ✅ resolvable despite the palindromic pair at the site |

The fourth row shows a strand flip selecting a different dosage, which only occurs in multiallelic matching. In the corresponding biallelic case, the complemented effect allele is absent. At `REF=A ALTS=[C,G]`, both readings select real alleles. The score does not establish its strand, so the matcher excludes it.

`--keep_ambiguous` deterministically selects the forward reading. A complemented reading without an other allele has no flipped match type. Therefore, `max()` sees one real match, not two. Exposure is small and decreasing because the PGS Catalog now requires `other_allele`. `REF_NO_OA` and `ALT_NO_OA` remain legacy paths.

**Exclude match collisions instead of choosing between them.** Two score variants in one accession can select the same allele of one target variant. Their weights then compete for one dosage. `any_value` could choose either weight and change the score, so the matcher excludes both score variants.

Best-match selection solves a different problem: several candidate matches for one score variant. `max(matches, key=match_priority)` selects among ALT candidates within one UDF call. A `ROW_NUMBER() OVER (PARTITION BY accession, row_nr, zarr_group ...)` window selects among target rows. v2 uses `_label_best_match` and `_label_duplicate_best_match` for this purpose and labels losing candidates `not_best`.

A collision instead involves several score variants selecting one target allele. Each row is already the best match for its own score variant. v2 handles this separately with `_label_duplicate_id`. Its docstring states: "refalt > altref > ... prioritisation doesn't fix this problem because row_nr is different" (`_match/label.py:178`).

v2 checks `(accession, ID)` over exploded biallelic variants. Multiallelic matching prevents reusing that key unchanged. v3 uses `(accession, target_row_nr, effect_allele_idx)`:

| two score variants, one accession, one position | v2 `(accession, ID)` | v3 `(accession, target_row_nr, effect_allele_idx)` |
|---|---|---|
| `A/C` and `A/G` → both index 0 | ❌ excluded | ❌ excluded — one dosage, two weights |
| `C/A` and `G/A` → index 1 and index 2 | ❌ excluded | ✅ both scored — two real alleles of one variant |
| `XY` and `X` canonicalised together → both index 0 | ❌ excluded | ❌ excluded — one dosage, two weights, reached by canonicalisation ([ADR 0008](0008-canonical-chromosome-names.md)) |

Multiallelic support creates the middle case. An exploded representation cannot retain it as two alleles of one site, because it gives them two IDs. Two score variants that match different alleles are two variants found, not one locus counted twice. Thus, matching each ALT does not inflate the `min_overlap` numerator.

Collisions can arise through these mechanisms:

- **Splitting multiallelics before the GWAS.** `bcftools norm -m -any` splits `A>C,G` into two records that both retain `REF=A`. Association results then contain two rows at one position. If both report the reference as their effect allele, both select index 0. This standard pipeline is the most likely source, not a malformed file.
- **Chromosome canonicalisation.** Converting `XY` to `X` places two score variants at one canonical position ([ADR 0008](0008-canonical-chromosome-names.md)).
- **Genome build remapping.** Two score variants can map to one harmonised position. v2's `_label_duplicate_id` docstring records `PGS002244` rows 1194115 and 1194132 both mapping to `3:50924580:C:A`.
- **Duplicates that the Catalog does not flag.** `check_duplicates` uses `chr:pos:effect_allele:other_allele` per accession for `is_duplicated` (`pgscatalog.core/lib/_normalise.py:89-90`). Rows that differ only in their other allele remain distinct by design.

Two score variants matching different ALTs at one multiallelic site are not a collision, as table row 2 shows. Nor is one score variant with several candidates. Best-match selection resolves those candidates before a result row exists.

**Do not add a `--keep_first_match` equivalent.** v2 provides one (`match_cli.py:246-252`). It could be added later if parity requires it. However, choosing the lowest `row_nr` is still arbitrary. v2's default also excludes every colliding candidate.

**Record a reason for every exclusion.** `EXCLUDED` covers four independent decisions. Without a reason, readers must infer the decision from `is_ambiguous` and `is_multiallelic`. Those flags describe the variant, not necessarily the exclusion.

For example, `is_multiallelic = true` describes the target on a matched row. On an excluded row, it can identify the cause. Add nullable `exclusion_reason` alongside those flags. It is NULL unless `match_summary = 'excluded'`.

| `exclusion_reason` | meaning | set by | suppressed by |
|---|---|---|---|
| `multiallelic` | target carries more than one alternative allele | `match_variant`, before any allele comparison | `--keep_multiallelic` |
| `ambiguous` | a complemented reading of the score variant also matched | `match_variant`, per score variant | `--keep_ambiguous` |
| `strand_flip` | the only reading that matched was a complemented one ([ADR 0011](0011-strand-flips-are-reverse-complements-and-optional.md)) | `match_variant`, per score variant | *not* passing `--ignore_strand_flips` |
| `collision` | another score variant of this accession resolved to the same target allele | the post-match pass, across rows | — |

`strand_flip` is the only exclusion enabled by request rather than by default. Its suppression rule therefore differs from the other flags.

The table defines precedence. `multiallelic` precedes allele comparisons and takes precedence over every later reason. `ambiguous` precedes selection of the winning match and takes precedence over `strand_flip`. `collision` comes last because it requires comparisons across score variants.

**This ADR rejects separate `match_type_enum` exclusion values.** The proposed values were `EXCLUDED_AMBIGUOUS`, `EXCLUDED_MULTIALLELIC` and `EXCLUDED_COLLISION`. `match_type` derives from `MatchPriority.name`. `MatchPriority` orders both `max()` and the ranking window. Extra exclusion members would require meaningless relative priorities in the type that controls match selection.

The summary groups by `exclusion_reason` and the existing flags. Instead of only reporting "217 variants excluded", it reports four counts with corresponding actions:

- Enable `--keep_multiallelic`.
- Enable `--keep_ambiguous`.
- Remove `--ignore_strand_flips`.
- Inspect the scoring file.

**bgen rejects multiallelic records instead of skipping them.** Format support is:

| input | biallelic | multiallelic | `*` / symbolic alleles | notes |
|---|---|---|---|---|
| VCF / BCF | ✅ | ✅ | ✅ kept positionally, never matched | no `--max-alleles 2` in the bcftools filter |
| bgen | ✅ | ❌ hard failure at LOAD | ❌ non-SNP records skipped | see below |
| plink (future reader) | ✅ | ❔ undecided | ❔ undecided | pvar comma-separated alts explode naturally |

`bgen_to_hardcall` assumes biallelic data through `probs[haploid, :2]` and a three-column unphased layout (`lib/cache/_bgen.py`). A triallelic, unphased diploid call has six probabilities in the specification's colex ordering. Support would need a general argmax and a decoder from index to allele pair. This is the largest part of the proposed work, for data that mostly does not exist.

The expected bgen inputs are population-scale biobank data from standard imputation pipelines. Those panels split multiallelics upstream. LOAD therefore rejects a multiallelic bgen record at a scoring-file position and names the variant. It must not silently score only a subset.

**This ADR rejects silently skipping these records.** Once VCF supports them, identical sites would behave differently by format without an explanation in the match log. **This ADR also rejects reading bgen alleles with undecodable genotypes.** That would require an "assessed but not decodable" Cache representation with no other use.

**`--keep_multiallelic` remains `false` by default**, as in v2 (`match_cli.py:238-241`). The default still improves reporting. Sites previously absent without explanation now receive `EXCLUDED` and `exclusion_reason = 'multiallelic'`. The old report said "not found" even though the site existed and policy discarded it.

A later release can change the default with one line, supported by release evidence and a changelog entry. This release must not combine that change with changes to Cache contents.

## Consequences

- **The excluded/unmatched distinction differs from v2.** v2 calculates ambiguity from `(ref, alt)` before checking for a match. It labels score variant `G/C` at palindromic target `A/T` as `EXCLUDED`, although neither allele matches. v3 labels it `NO_MATCH`. Scores stay unchanged because `match_rate` counts only `matched`, but log comparisons expose the difference.
- **Allele index capacity is sufficient.** Genotypes use `uint8`, with 254 and 255 reserved (`lib/constants.py`). A variant can therefore have 253 ALTs before the encoding fails.
- **Earlier caches contain a filtered variant set.** They represent multiallelic positions as absent padding. As with [ADR 0008](0008-canonical-chromosome-names.md), discard alpha caches instead of migrating them.
