# Multiallelic variants are matched per alternative allele; bgen stays biallelic

## Status

Accepted. Its definition of a complemented reading is refined by [ADR 0011](0011-a-strand-flip-is-a-reverse-complement-and-indels-decline-it.md).

## Context

`pgscatalog.match` (v2) supports only the first alternative allele because it has to hand plink2 a biallelic pvar: `annotate_multiallelic` splits `ALT` on commas and explodes the frame (`_match/preprocess.py:42-59`), then matching looks for the biallelic equivalent of each scoring variant. The representation is lossy, and the natural worry carried over into v3 is what a strand flip does to the dosage of an exploded allele.

In v3 that worry is unfounded, and for a reason worth writing down. `match_variant` returns an `effect_allele_idx`, and `calculate_effect_allele_dosage` counts genotype calls equal to that index (`lib/score/_dosage.py:317`). Allele *sequences* are compared once, at match time; dosage is index arithmetic afterwards. Complementing `A`→`T` does not move an allele's position in the record, so `REFALT_FLIP` and `ALTREF_FLIP` yield the same index as their unflipped counterparts. A flip changes *which* allele a score matched, never how a matched allele is counted. v3 needs no explode step and should not grow one.

What blocks multiallelic support is therefore not the dosage path but the cache boundary, plus three defects behind it:

- **Nothing multiallelic reaches the cache.** VCF ingest passes `--max-alleles 2` (`lib/cache/_vcf.py:262-265`) and bgen skips (`lib/cache/_bgen.py:283-285`). `--keep_multiallelic` is unreachable code. Excluded sites are padded with `ref=None`, filtered out by `update_match_table`'s `WHERE target.ref IS NOT NULL`, and surface in the match log as ordinary unmatched variants — the user is told a variant was not found when it was found and deliberately discarded.
- **Only `alts[0]` is examined.** `get_variant_match_priority` breaks after the first alternative allele when multiallelic (`lib/score/_matchvariants.py:355-357`), but the ambiguous branch `continue`s past that break (line 338). `REF=A ALTS=[C,T]` against a score variant `A/T` reports *unmatched*; the same site written `ALTS=[T,C]` reports *excluded*. The answer depends on ALT column order, which carries no meaning.
- **Non-ACGT alleles raise inside a UDF.** `is_biallelic_variant_ambiguous` raises `ValueError` on anything outside `ACGT` (line 371). This was believed unreachable. It is not: `--exclude-types other` excludes nothing (verified against bcftools — a `G,<NON_REF>` record passes it untouched), so `--max-alleles 2` is the only filter in play, and a biallelic `A/*` record passes it. `match_variant(..., ref='A', alts=['*'])` raises today, on `main`.

## Decision

**Multiallelic records are ingested, and every alternative allele is a match candidate.** The per-alt `break` goes; a score variant is matched against each alternative allele and the best priority wins. This needs no tie-break policy, because at most one alternative allele can match a given score variant: with an other allele, a forward match on `alt_i` and a flipped match on `alt_j` both require `alt_i = complement(ref) = alt_j`; without one, the effect allele equals exactly one distinct alternative allele.

**Alleles are positional, so unmatchable alleles are kept and skipped, never filtered.** Genotype calls name alleles by index, so dropping `*` from `[C,*,G]` renumbers `G` from 3 to 2 and every dosage for it becomes the wrong allele. The spanning-deletion `*` and symbolic alleles like `<NON_REF>` stay in the list, are skipped at match time without an ambiguity check, and yield `NO_MATCH`. `REF=A ALTS=[*,G]` against `A/G` matches `G` at index 2 and scores normally. The `ValueError` survives for genuinely corrupt allele strings only.

**Ambiguity is redefined as a property of the pairing, not of a pair of alleles.** A score variant is ambiguous when *a complemented reading of it also matches an allele the target carries*. This is a strict generalisation of the palindrome check, so `is_biallelic_variant_ambiguous` and the `match_ambiguous` pre-check are deleted and ambiguity is computed from whether a second reading matched. What "a complemented reading" means is settled by [ADR 0011](0011-a-strand-flip-is-a-reverse-complement-and-indels-decline-it.md), written after this one: it is a *reverse* complement, and an indel gets no complemented reading at all. Every row of the table below is single-base, where the two definitions agree, so none of it moves:

| target | score variant | forward | complemented | verdict |
|---|---|---|---|---|
| `A / T` | `A/T` | idx 0 | idx 1 | ✅ ambiguous — palindrome, as before |
| `A / C` | `A/C` | idx 0 | — | ✅ resolvable |
| `A / C` | `C` (no other allele) | idx 1 | — | ✅ resolvable |
| `A / C,G` | `C` (no other allele) | idx 1 | idx 2 | ✅ ambiguous — only reachable multiallelically |
| `A / T,C` | `C/A` | idx 2 | — | ✅ resolvable despite the palindromic pair at the site |

The fourth row is the one case where a strand flip really does change which dosage gets counted, and it is new. `classify_match` attempts no complemented reading when the other allele is NULL (`lib/score/_matchvariants.py:304-309`), which is harmless biallelically — the complemented effect allele is not at the site, so there is nothing to confuse it with. At `REF=A ALTS=[C,G]` both readings land on real alleles and the forward one wins silently. Under the new definition it is excluded instead. The cost is small and shrinking: the PGS Catalog now enforces `other_allele`, so `REF_NO_OA`/`ALT_NO_OA` are legacy paths.

**Match collisions are excluded, not resolved.** Two score variants of one accession can resolve to the same allele of the same target variant, and then two weights compete for one dosage. `create_wide_weights_table`'s PIVOT settles it with `coalesce(any_value(effect_weight), 0)` (`lib/score/_weight_matrix.py:116`) — a coin flip that changes the score, beside a comment asserting the very uniqueness that has been violated. Colliding score variants are excluded instead.

Best-match selection cannot fix this and it is worth being explicit about why, because the two look alike. Best-match runs *many candidate matches → one score variant*: `max(matches, key=match_priority)` (`lib/score/_matchvariants.py:208`) collapses the per-alt candidates inside a single UDF call, and the `ROW_NUMBER() OVER (PARTITION BY accession, row_nr, zarr_group ...)` window (lines 107-114) collapses candidates across target rows. v2 does the same thing in `_label_best_match` and `_label_duplicate_best_match`, labelling the losers `not_best`. A collision runs the other way — *many score variants → one target allele* — so nothing is competing to be best; each row already **is** the best match for its own score variant. v2's `_label_duplicate_id` handles this case separately, and its docstring says exactly this: "refalt > altref > ... prioritisation doesn't fix this problem because row_nr is different" (`_match/label.py:178`).

v2's key for that check is `(accession, ID)` over exploded biallelic variants, and it **cannot be ported unchanged** — multiallelic matching is precisely what breaks it. The key here is `(accession, target_row_nr, effect_allele_idx)`:

| two score variants, one accession, one position | v2 `(accession, ID)` | v3 `(accession, target_row_nr, effect_allele_idx)` |
|---|---|---|
| `A/C` and `A/G` → both index 0 | ❌ excluded | ❌ excluded — one dosage, two weights |
| `C/A` and `G/A` → index 1 and index 2 | ❌ excluded | ✅ both scored — two real alleles of one variant |
| `XY` and `X` canonicalised together → both index 0 | ❌ excluded | ❌ excluded — one dosage, two weights, reached by canonicalisation ([ADR 0008](0008-canonical-chromosome-names.md)) |

The middle row is the state multiallelic support creates and an exploded representation cannot reach, because there two alleles of one site are two IDs.

What produces a collision, so a reader can recognise one:

- **Split multiallelics upstream of the GWAS.** `bcftools norm -m -any` splits `A>C,G` into two records that both keep `REF=A`. Association results computed on the split panel are two rows at one position, and if the effect allele is reported as the reference in both, both resolve to index 0. This is the mechanism most likely to be met in the wild, and it is a standard pipeline, not a malformed file.
- **Chromosome canonicalisation.** `XY` collapsing into `X` ([ADR 0008](0008-canonical-chromosome-names.md)) puts two score variants at one canonical position.
- **Genome build remapping.** Two score variants remapped onto one harmonised position — v2's `_label_duplicate_id` docstring records a real instance, `PGS002244` rows 1194115 and 1194132 both landing on `3:50924580:C:A`.
- **Duplicate rows the Catalog does not flag.** `check_duplicates` keys `is_duplicated` on `chr:pos:effect_allele:other_allele` per accession (`pgscatalog.core/lib/_normalise.py:89-90`), so two rows differing only in their other allele pass through as distinct by design.

What is *not* a collision: two score variants matching different alternative alleles of one multiallelic variant (row 2 of the table) — two variants found, not one locus counted twice — and a single score variant producing several candidate matches, which best-match selection collapses before a row ever exists.

**No `--keep_first_match` equivalent is added.** v2 offers one (`match_cli.py:247-252`) and it can be added later if parity turns out to matter, but keeping the lowest `row_nr` is the same coin flip with a nicer name, and v2's own default is to drop every colliding candidate.

**Every exclusion names its reason.** `EXCLUDED` now covers three unrelated decisions, and today a reader has to infer which from `is_ambiguous` and `is_multiallelic` — flags that describe the *variant*, not the *decision*. They are not interchangeable: `is_multiallelic = true` on a matched row is a fact about the target, on an excluded row it is the cause. A nullable `exclusion_reason` column is added alongside them, NULL unless `match_summary = 'excluded'`:

| `exclusion_reason` | meaning | set by | suppressed by |
|---|---|---|---|
| `multiallelic` | target carries more than one alternative allele | `match_variant`, before any allele comparison | `--keep_multiallelic` |
| `ambiguous` | a complemented reading of the score variant also matched | `match_variant`, per alternative allele | `--keep_ambiguous` |
| `collision` | another score variant of this accession resolved to the same target allele | the post-match pass, across rows | — |

`multiallelic` takes precedence over `ambiguous` where both apply, which the existing code order already gives for free — `match_variant` tests `is_multiallelic and not match_multiallelic` before matching any alleles (`lib/score/_matchvariants.py:190`). `collision` is necessarily last, being the only reason that is not a property of one score variant.

**Splitting `match_type_enum` into `EXCLUDED_AMBIGUOUS` / `EXCLUDED_MULTIALLELIC` / `EXCLUDED_COLLISION` was rejected.** `match_type` is derived from `MatchPriority.name`, and `MatchPriority` is an ordering used by `max()` and by the ranking window. Three exclusion members would need three priority values whose relative order means nothing, encoding a false ranking into the type that drives match selection.

The summary log groups on `exclusion_reason` alongside the existing flags, so "217 variants excluded" becomes three counts that each name a fix — raise `--keep_multiallelic`, raise `--keep_ambiguous`, or go and look at the scoring file.

**bgen refuses multiallelic records rather than skipping them.** Format support is therefore:

| input | biallelic | multiallelic | `*` / symbolic alleles | notes |
|---|---|---|---|---|
| VCF / BCF | ✅ | ✅ | ✅ kept positionally, never matched | `--max-alleles 2` removed from the bcftools filter |
| bgen | ✅ | ❌ hard failure at LOAD | ❌ non-SNP records skipped, as now | see below |
| plink (future reader) | ✅ | ❔ undecided | ❔ undecided | pvar comma-separated alts explode naturally |

`bgen_to_hardcall` is hard-wired biallelic — `probs[haploid, :2]` and a three-column unphased layout (`lib/cache/_bgen.py:164-170`). A triallelic unphased diploid call carries six probabilities in the spec's colex ordering, so support means a general argmax plus an index→allele-pair decoder: the largest single piece of work in this feature, for data that mostly does not exist. bgen input here is population-scale biobank data from standard imputation pipelines, and those panels split multiallelics upstream.

The `NotImplementedError` at `lib/cache/_bgen.py:286-287` therefore fires on the skip path too, and `skip_multiallelic` is deleted. A bgen holding a multiallelic record at a scoring-file position refuses to load, naming the variant, rather than quietly scoring a subset. **Silently skipping was rejected**: once VCF matches multiallelics, identical sites would behave differently by input format with nothing in the match log saying which happened. **Ingesting bgen alleles with undecodable genotypes was also rejected** — it needs a cache representation for "assessed but not decodable" that nothing else wants.

**`--keep_multiallelic` keeps its default of `false`**, matching v2 (`match_cli.py:233-236`). The default path still improves: those sites go from silently absent to `EXCLUDED` with `is_multiallelic = true` in the match log. Flipping the default is a one-line change to make later, on release evidence, with a changelog entry — not something to slip into a release that also changes what the cache holds.

## Consequences

- **The fix first proposed for a duplicate canonical position becomes wrong.** It suggested a `(accession, chr_name, chr_position)` uniqueness check on `score_variant_table`. That key rejects valid scoring files once multiallelics land, and it sits on the wrong table — a collision is a post-match property, invisible while the scoring file is being read. What that fix had to achieve is unchanged; the key it proposed is superseded by the collision key above.
- **0027's second wrong is narrower than it looks.** It reports the `min_overlap` numerator counting one locus twice. Two score variants matching two *different* alleles of one multiallelic variant are two variants found, not one locus double-counted, so the inflation is real only for same-index collisions — which excluding collisions already removes.
- **The excluded/unmatched split in every match log shifts.** Ambiguity is currently computed from `(ref, alt)` before anyone checks whether the score variant matched, so a score variant of `G/C` at a palindromic `A/T` site — matching nothing — is logged `EXCLUDED` today and becomes `NO_MATCH` under the new definition. No score changes: `match_rate` sums only `matched` (`lib/score/_matchlog.py:137`). Anyone diffing v3 logs against v2 will see it.
- **The `A/*` crash is fixed as a side effect**, not as separate work. It is a pre-existing bug on `main`, reachable without any multiallelic record present.
- **`ref_first` is deleted** (`lib/cache/_bgen.py:254, 297-299`). It swaps which allele the cache labels REF, but `bgen_to_hardcall` returns argmax indices in the bgen's own allele order, so `ref_first=False` makes a `REFALT` match with `effect_allele_idx=0` count the alt. Nothing passes it and no test covers it. It is the same failure this ADR rules out for strand, arriving through allele ordering instead.
- **Allele index headroom is not a constraint.** Genotypes are `uint8` with 254 and 255 reserved (`lib/constants.py:9-18`), so a variant can carry 253 alternative alleles before the encoding breaks.
- **Caches written before this decision hold a filtered variant set** — multiallelic positions padded as absent. As with [ADR 0008](0008-canonical-chromosome-names.md), alpha caches are discarded rather than migrated.
