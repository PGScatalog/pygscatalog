# Multiallelic variants are matched per alternative allele; bgen stays biallelic

## Status

Accepted. Covers how a score variant is matched against a target variant carrying
more than one alternative allele. What a complemented reading of a score variant is,
and when one is attempted, is [ADR 0011](0011-strand-flips-are-reverse-complements-and-optional.md).

## Context

`pgscatalog.match` (v2) supports only the first alternative allele because it has to
hand plink2 a biallelic pvar: `annotate_multiallelic` splits `ALT` on commas and
explodes the frame (`_match/preprocess.py:42-59`), then matching looks for the
biallelic equivalent of each scoring variant. The representation is lossy, and the
natural worry carried over into v3 is what a strand flip does to the dosage of an
exploded allele.

In v3 that worry is unfounded, and for a reason worth writing down. `match_variant`
returns an `effect_allele_idx`, and `calculate_effect_allele_dosage` counts genotype
calls equal to that index (`lib/score/_dosage.py`). Allele *sequences* are compared
once, at match time; dosage is index arithmetic afterwards. Complementing `A`→`T`
does not move an allele's position in the record, so `REFALT_FLIP` and `ALTREF_FLIP`
yield the same index as their unflipped counterparts. A flip changes *which* allele a
score matched, never how a matched allele is counted. v3 needs no explode step and
should not grow one.

What multiallelic support costs instead is at the cache boundary: the bcftools filter
has to stop dropping these records, bgen has to say what it does with one, and the
matcher has to cope with alleles that are not sequences at all. `--exclude-types
other` excludes nothing (verified against bcftools — a `G,<NON_REF>` record passes it
untouched), so a record can reach the matcher carrying `*` or `<NON_REF>` whether or
not it is multiallelic.

## Decision

**Multiallelic records are ingested, and every alternative allele is a match
candidate.** A score variant is matched against each alternative allele and the best
priority wins. This needs no tie-break policy, because at most one alternative allele
can match a given score variant: with an other allele, a forward match on `alt_i` and
a complemented match on `alt_j` both require `alt_i = reverse_complement(ref) =
alt_j`; without one, the effect allele equals exactly one distinct alternative allele.

**Alleles are positional, so unmatchable alleles are kept and skipped, never
filtered.** Genotype calls name alleles by index, so dropping `*` from `[C,*,G]`
renumbers `G` from 3 to 2 and every dosage for it becomes the wrong allele. The
spanning-deletion `*` and symbolic alleles like `<NON_REF>` stay in the list, are
skipped at match time without an ambiguity check, and yield `NO_MATCH`. `REF=A
ALTS=[*,G]` against `A/G` matches `G` at index 2 and scores normally. Anything else
that is not a nucleotide sequence is a corrupt allele string and raises.

**Ambiguity is a property of the pairing, not of a pair of alleles.** A score variant
is ambiguous when *a complemented reading of it also matches an allele the target
carries*. This is a strict generalisation of the palindrome check, so there is no
separate palindrome test: ambiguity is computed from whether a second reading
matched. Every row below is single-base, where a plain complement and a reverse
complement agree, so none of it depends on which is meant:

| target | score variant | forward | complemented | verdict |
|---|---|---|---|---|
| `A / T` | `A/T` | idx 0 | idx 1 | ✅ ambiguous — palindrome |
| `A / C` | `A/C` | idx 0 | — | ✅ resolvable |
| `A / C` | `C` (no other allele) | idx 1 | — | ✅ resolvable |
| `A / C,G` | `C` (no other allele) | idx 1 | idx 2 | ✅ ambiguous — only reachable multiallelically |
| `A / T,C` | `C/A` | idx 2 | — | ✅ resolvable despite the palindromic pair at the site |

The fourth row is the one case where a strand flip really does change which dosage
gets counted, and it is only reachable multiallelically: biallelically the
complemented effect allele is not at the site, so there is nothing to confuse it
with. At `REF=A ALTS=[C,G]` both readings land on real alleles, so the score variant
does not say which strand it was written on and it is excluded. `--keep_ambiguous`
resolves it to the forward reading, deterministically — a no-OA complemented reading
has no flipped match type to be promoted to, so `max()` sees one real match, not two.
The cost is small and shrinking: the PGS Catalog now enforces `other_allele`, so
`REF_NO_OA`/`ALT_NO_OA` are legacy paths.

**Match collisions are excluded, not resolved.** Two score variants of one accession
can resolve to the same allele of the same target variant, and then two weights
compete for one dosage. Picking one with `any_value` is a coin flip that changes the
score, so colliding score variants are excluded instead.

Best-match selection cannot fix this and it is worth being explicit about why,
because the two look alike. Best-match runs *many candidate matches → one score
variant*: `max(matches, key=match_priority)` collapses the per-alt candidates inside a
single UDF call, and a `ROW_NUMBER() OVER (PARTITION BY accession, row_nr,
zarr_group ...)` window collapses candidates across target rows. v2 does the same
thing in `_label_best_match` and `_label_duplicate_best_match`, labelling the losers
`not_best`. A collision runs the other way — *many score variants → one target
allele* — so nothing is competing to be best; each row already **is** the best match
for its own score variant. v2's `_label_duplicate_id` handles this case separately,
and its docstring says exactly this: "refalt > altref > ... prioritisation doesn't fix
this problem because row_nr is different" (`_match/label.py:178`).

v2's key for that check is `(accession, ID)` over exploded biallelic variants, and it
**cannot be ported unchanged** — multiallelic matching is precisely what breaks it.
The key here is `(accession, target_row_nr, effect_allele_idx)`:

| two score variants, one accession, one position | v2 `(accession, ID)` | v3 `(accession, target_row_nr, effect_allele_idx)` |
|---|---|---|
| `A/C` and `A/G` → both index 0 | ❌ excluded | ❌ excluded — one dosage, two weights |
| `C/A` and `G/A` → index 1 and index 2 | ❌ excluded | ✅ both scored — two real alleles of one variant |
| `XY` and `X` canonicalised together → both index 0 | ❌ excluded | ❌ excluded — one dosage, two weights, reached by canonicalisation ([ADR 0008](0008-canonical-chromosome-names.md)) |

The middle row is the state multiallelic support creates and an exploded
representation cannot reach, because there two alleles of one site are two IDs. It is
also why matching per alternative allele does not inflate the `min_overlap`
numerator: two score variants matching two *different* alleles of one multiallelic
variant are two variants found, not one locus counted twice.

What produces a collision, so a reader can recognise one:

- **Split multiallelics upstream of the GWAS.** `bcftools norm -m -any` splits
  `A>C,G` into two records that both keep `REF=A`. Association results computed on
  the split panel are two rows at one position, and if the effect allele is reported
  as the reference in both, both resolve to index 0. This is the mechanism most
  likely to be met in the wild, and it is a standard pipeline, not a malformed file.
- **Chromosome canonicalisation.** `XY` collapsing into `X`
  ([ADR 0008](0008-canonical-chromosome-names.md)) puts two score variants at one
  canonical position.
- **Genome build remapping.** Two score variants remapped onto one harmonised
  position — v2's `_label_duplicate_id` docstring records a real instance,
  `PGS002244` rows 1194115 and 1194132 both landing on `3:50924580:C:A`.
- **Duplicate rows the Catalog does not flag.** `check_duplicates` keys
  `is_duplicated` on `chr:pos:effect_allele:other_allele` per accession
  (`pgscatalog.core/lib/_normalise.py:89-90`), so two rows differing only in their
  other allele pass through as distinct by design.

What is *not* a collision: two score variants matching different alternative alleles
of one multiallelic variant (row 2 of the table) — two variants found, not one locus
counted twice — and a single score variant producing several candidate matches, which
best-match selection collapses before a row ever exists.

**No `--keep_first_match` equivalent is added.** v2 offers one
(`match_cli.py:246-252`) and it can be added later if parity turns out to matter, but
keeping the lowest `row_nr` is the same coin flip with a nicer name, and v2's own
default is to drop every colliding candidate.

**Every exclusion names its reason.** `EXCLUDED` covers four unrelated decisions, and
without a reason column a reader has to infer which from `is_ambiguous` and
`is_multiallelic` — flags that describe the *variant*, not the *decision*. They are
not interchangeable: `is_multiallelic = true` on a matched row is a fact about the
target, on an excluded row it is the cause. A nullable `exclusion_reason` column sits
alongside them, NULL unless `match_summary = 'excluded'`:

| `exclusion_reason` | meaning | set by | suppressed by |
|---|---|---|---|
| `multiallelic` | target carries more than one alternative allele | `match_variant`, before any allele comparison | `--keep_multiallelic` |
| `ambiguous` | a complemented reading of the score variant also matched | `match_variant`, per score variant | `--keep_ambiguous` |
| `strand_flip` | the only reading that matched was a complemented one ([ADR 0011](0011-strand-flips-are-reverse-complements-and-optional.md)) | `match_variant`, per score variant | *not* passing `--ignore_strand_flips` |
| `collision` | another score variant of this accession resolved to the same target allele | the post-match pass, across rows | — |

`strand_flip` is the only row that exists by request rather than by default, which is
why its "suppressed by" reads backwards from the others.

Precedence runs down the table. `multiallelic` is tested before any allele is
compared, so it wins over everything below it. `ambiguous` is decided before a
winning match is known, so it wins over `strand_flip`. `collision` is necessarily
last, being the only reason that is not a property of one score variant.

**Splitting `match_type_enum` into `EXCLUDED_AMBIGUOUS` / `EXCLUDED_MULTIALLELIC` /
`EXCLUDED_COLLISION` was rejected.** `match_type` is derived from
`MatchPriority.name`, and `MatchPriority` is an ordering used by `max()` and by the
ranking window. Extra exclusion members would need priority values whose relative
order means nothing, encoding a false ranking into the type that drives match
selection.

The summary log groups on `exclusion_reason` alongside the existing flags, so "217
variants excluded" becomes four counts that each name a fix — raise
`--keep_multiallelic`, raise `--keep_ambiguous`, drop `--ignore_strand_flips`, or go
and look at the scoring file.

**bgen refuses multiallelic records rather than skipping them.** Format support is
therefore:

| input | biallelic | multiallelic | `*` / symbolic alleles | notes |
|---|---|---|---|---|
| VCF / BCF | ✅ | ✅ | ✅ kept positionally, never matched | no `--max-alleles 2` in the bcftools filter |
| bgen | ✅ | ❌ hard failure at LOAD | ❌ non-SNP records skipped | see below |
| plink (future reader) | ✅ | ❔ undecided | ❔ undecided | pvar comma-separated alts explode naturally |

`bgen_to_hardcall` is hard-wired biallelic — `probs[haploid, :2]` and a three-column
unphased layout (`lib/cache/_bgen.py`). A triallelic unphased diploid call carries six
probabilities in the spec's colex ordering, so support means a general argmax plus an
index→allele-pair decoder: the largest single piece of work in this feature, for data
that mostly does not exist. bgen input here is population-scale biobank data from
standard imputation pipelines, and those panels split multiallelics upstream. A bgen
holding a multiallelic record at a scoring-file position therefore refuses to load,
naming the variant, rather than quietly scoring a subset.

**Silently skipping was rejected**: once VCF matches multiallelics, identical sites
would behave differently by input format with nothing in the match log saying which
happened. **Ingesting bgen alleles with undecodable genotypes was also rejected** — it
needs a cache representation for "assessed but not decodable" that nothing else
wants.

**`--keep_multiallelic` keeps its default of `false`**, matching v2
(`match_cli.py:238-241`). The default path still improves: those sites go from
silently absent to `EXCLUDED` with `exclusion_reason = 'multiallelic'` in the match
log, where before the user was told a variant was not found when it was found and
deliberately discarded. Flipping the default is a one-line change to make later, on
release evidence, with a changelog entry — not something to slip into a release that
also changes what the cache holds.

## Consequences

- **The excluded/unmatched split differs from v2.** v2 computes ambiguity from
  `(ref, alt)` before anyone checks whether the score variant matched, so a score
  variant of `G/C` at a palindromic `A/T` site — matching nothing — is `EXCLUDED`
  there and `NO_MATCH` here. No score changes: `match_rate` sums only `matched`.
  Anyone diffing v3 logs against v2 will see it.
- **Allele index headroom is not a constraint.** Genotypes are `uint8` with 254 and
  255 reserved (`lib/constants.py`), so a variant can carry 253 alternative alleles
  before the encoding breaks.
- **Caches written before this decision hold a filtered variant set** — multiallelic
  positions padded as absent. As with [ADR 0008](0008-canonical-chromosome-names.md),
  alpha caches are discarded rather than migrated.
