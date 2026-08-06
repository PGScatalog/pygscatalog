# A strand flip is a reverse complement, and an indel declines it

## Status

Accepted. Corrects the ambiguity section of [ADR 0010](0010-multiallelic-matched-per-alt-bgen-stays-biallelic.md), which defines ambiguity in terms of "a complemented reading" without saying what one is.

## Context

`complement` translated `ATCG` in place and stopped there, and `classify_match` applied it at any allele length. The strand-flipped reading of a sequence is its **reverse** complement: a sequence is read 5'→3' on both strands, so flipping strand reverses the order of the bases as well as complementing each one.

```
'A'    complement 'T'     reverse complement 'T'    <- agree
'AT'   complement 'TA'    reverse complement 'AT'   <- differ
'GC'   complement 'CG'    reverse complement 'GC'
'ATG'  complement 'TAC'   reverse complement 'CAT'
```

For a single base the reversal is invisible, which is why nothing noticed. For anything longer the two disagree, so every `_FLIP` match type matched the sequences a strand flip does not produce and missed the ones it does: a score variant `AC/CC` was reported `REFALT_FLIP` against a target `TG/GG`, which is not what the other strand of `AC/CC` carries, and reported unmatched against `GT/GG`, which is.

This predates the per-alternative-allele matching of [ADR 0010](0010-multiallelic-matched-per-alt-bgen-stays-biallelic.md), which widened the blast radius rather than creating it. The old `is_biallelic_variant_ambiguous` returned `False` for anything longer than one base, so a multi-base pairing could never be called ambiguous. ADR 0010 redefined ambiguity as "a complemented reading of the score variant also matched an allele the target carries", so the same wrong complement began deciding exclusions as well as matches.

`pgscatalog.match` (v2) has the identical defect. `complement_valid_alleles` (`_match/preprocess.py:18`) does placeholder character substitution with no reversal, guarded only by `^[ACGT]+$` and applied to indels as well as SNPs; `get_all_matches` attempts all four flip strategies for every variant with no length or SNP restriction (`_match/match.py:44-67`); and `_label_biallelic_ambiguous` (`_match/label.py:223`) uses the same plain complement with no length guard, so `AT`/`TA` is flagged ambiguous there too.

plink2 contributes no strand logic to compare against. v2's flipped allele becomes `matched_effect_allele` and is renamed `effect_allele` in the plink2 scorefile (`_match/plink.py:96`), so plink2 matches a literal string against the pvar and has no opinion about strand. There is no external behaviour to stay compatible with: the decision is entirely ours.

Indels are the reason this needed a decision rather than a one-line fix. A VCF indel allele carries a shared anchor base at the start (`REF=A ALT=AT`), and that anchor is *positional* — it stays leftmost on whichever strand the record is read.

## Decision

**A complemented reading is a reverse complement.** `reverse_complement` complements each base and reverses the sequence. Single-base alleles are unchanged from before, so SNPs — nearly all of the exposure — do not move.

**No complemented reading is attempted where the two alleles differ in length.** An indel can still match forward (`REFALT`, `ALTREF`, `REF_NO_OA`, `ALT_NO_OA`); it simply never matches via a complemented reading, and is therefore never ambiguous either.

The rejected alternative is reverse-complementing indels anyway. Because the anchor base is positional, the strand-flipped record sits at a different position with different alleles, and no string operation on the two alleles in hand recovers it. Reverse-complementing an indel produces something plausible-looking and wrong, which is exactly the failure mode this ADR exists to remove. Declining is the honest answer: an indel written on the other strand is a record this matcher cannot recognise, and saying so is better than matching the wrong one.

**The length comparison is made between the target's reference allele and the one alternative allele being read**, and it gates the whole complemented reading rather than one branch of it. "The two alleles are the same length" is a property of the *pairing* being read, never of one allele alone: `A` says nothing about whether the site is an indel, and a score variant with no other allele has no second allele of its own to compare with. The target pair is the one that is always available — biallelic or multiallelic, with an other allele or without — and it is per alternative allele, so a multiallelic variant carrying both a SNP and an indel alternative gets the right answer for each.

It is also sufficient for the readings that *do* have a score-side pair, because a reverse complement preserves length: a complemented match against an equal-length target pair forces the score alleles to be equal-length too. Checking the score side as well would add a condition that can never fire.

**The SQL `complement` macro becomes `reverse_complement` and loses its `upper()`.** `score_log_table` derives `matched_effect_allele` with the macro while matching is decided by the Python function, so the two have to be the same function or the log names a different allele from the match that produced it. The `upper()` was never load-bearing: every allele comparison in the matcher is case-sensitive and `check_nucleotide_allele` refuses a target allele that is not uppercase `ACGT`, so no row carrying a lowercase allele can reach a flipped match type. What `upper()` could do was print a sequence the matcher never compared, so it goes.

**`matched_effect_allele` for `ALTREF_FLIP` is corrected to read the effect allele.** It read `complement(other_allele)`, which resolves to the target's *reference* allele — the allele that pairing did not match. Both flipped match types resolve the score's effect allele to `reverse_complement(effect_allele)`; they differ in *which* target allele that turns out to be, not in which score allele was read. v2 gets this right (`effect_allele_FLIP` for both `refalt_flip` and `altref_flip`), so this is v3 catching up rather than diverging. It was found while making the two implementations of the same function agree, and fixed here because leaving it would have meant testing an assertion known to be wrong.

## Consequences

- **v3 knowingly diverges from v2 on multi-base alleles, in cases where v2 is wrong.** A multi-base pairing that v2 reports `REFALT_FLIP` may be unmatched here, and one v2 reports unmatched may be a flipped match here. Anyone diffing v3 match logs against v2 will see it, and the difference is the fix.
- **Multi-base `_FLIP` matches change in both directions.** This is not purely additive: MNV pairings that matched a plain complement stop matching, and the genuine strand flips they were standing in for start matching. `match_flipped` counts in the summary log move accordingly. Single-base alleles are byte-identical, so the practical exposure is small.
- **Ambiguity moves with it**, because ADR 0010 computes ambiguity from whether a complemented reading matched. Two changes fall out. A multi-base pairing is now ambiguous only when it is genuinely palindromic under reverse complement — `AT`/`TA` still is, since `AT` and `TA` are each their own reverse complement, but it is so for the right reason now. And an indel is never ambiguous at all, so an indel that matches forward is scored rather than excluded.
- **An indel's strand is now unresolvable rather than wrongly resolved.** Scoring files carrying indels written on the opposite strand from the target lose those variants to `NO_MATCH`, where before some of them were matched to the wrong allele. That costs match rate and buys correctness; there is no representation of a strand-flipped indel that the two alleles alone can reach.
- **`complement` is gone as a name, in Python and in SQL.** The name was the bug: it did what it said and what it said was not a strand flip. `complement_table` survives as what it is, the translation table `reverse_complement` is built from.
- **ADR 0010's ambiguity table stays correct.** Every row in it is single-base, so a plain complement and a reverse complement agree throughout; that ADR's reasoning does not depend on which one is meant, but the phrase "a complemented reading" now has a definition, and it is this one.
