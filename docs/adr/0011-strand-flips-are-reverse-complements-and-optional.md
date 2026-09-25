# A strand flip is a reverse complement, and matching one is optional

## Status

Accepted. Defines the complemented reading used to calculate ambiguity in [ADR 0010](0010-multiallelic-matched-per-alt-bgen-stays-biallelic.md). Adds a flag to reject matches that require that reading.

## Context

A scoring file and a target genome can report the same variant from opposite strands. The strand-flipped reading of a sequence is its **reverse complement**. Both strands are read 5'→3', so changing strand reverses the base order and complements each base.

```
'A'    complement 'T'     reverse complement 'T'    <- agree
'AT'   complement 'TA'    reverse complement 'AT'   <- differ
'GC'   complement 'CG'    reverse complement 'GC'
'ATG'  complement 'TAC'   reverse complement 'CAT'
```

Reversal has no visible effect on a single base. SNPs account for nearly all cases, so a plain complement can appear correct. Longer sequences differ. A plain complement matches sequences that a strand flip does not produce and misses sequences that it does produce.

For example, a plain complement labels score variant `AC/CC` as `REFALT_FLIP` against target `TG/GG`. That target does not represent the opposite strand of `AC/CC`. The correct target is `GT/GG`, which the plain complement leaves unmatched.

Indels require a separate decision. A VCF indel has a shared anchor base at the start, as in `REF=A ALT=AT`. This anchor is positional and stays leftmost on either strand.

The complemented reading serves two purposes. It determines `REFALT_FLIP` and `ALTREF_FLIP` matches. [ADR 0010](0010-multiallelic-matched-per-alt-bgen-stays-biallelic.md) also uses it to detect ambiguity. A score variant is ambiguous when its complemented reading also matches an allele in the target.

plink2 provides no strand logic for this comparison. v2 records the flipped allele as `matched_effect_allele`. It renames that field to `effect_allele` in the plink2 scorefile (`_match/plink.py:96`). plink2 then compares a literal string against the pvar, without interpreting strand. Thus, no external strand behaviour constrains this decision.

Some users know the strand of their target genotypes. They prefer to discard a variant rather than infer its strand. v2 provides `--ignore_strand_flips` for them (`match_cli.py:246`). v3 also needs a policy for these users.

## Decision

**A complemented reading is a reverse complement.** `reverse_complement` complements each base and reverses the sequence. The SQL macro that derives `matched_effect_allele` for the log uses the same function. The log therefore names the allele selected by the match.

**Do not attempt a complemented reading when the two alleles differ in length.** An indel can still match forward through `REFALT`, `ALTREF`, `REF_NO_OA` or `ALT_NO_OA`. It cannot match through a complemented reading and cannot be ambiguous.

This ADR rejects reverse-complementing indels. The positional anchor means that the opposite-strand record has a different position and different alleles. No string operation on the available allele pair can reconstruct that record. Reverse-complementing the indel can instead produce a plausible but incorrect result.

The matcher must report its limitation rather than select an incorrect allele. It cannot recognise an indel written on the opposite strand through those two allele strings alone.

**Compare the target REF length with the length of the ALT under consideration.** This comparison controls the whole complemented reading, not one branch. Equal length describes an allele pair, not an individual allele. An `A` alone does not establish whether the variant is an indel. A score variant without an other allele provides no score-side pair for comparison.

The target pair always exists, regardless of allele count or whether the score provides an other allele. The check runs separately for each ALT. Thus, a multiallelic target with SNP and indel alternatives receives the correct policy for each alternative.

A separate score-side length check would add no protection. Reverse complements preserve length. A complemented match to an equal-length target pair therefore requires equal-length score alleles too.

**`--ignore_strand_flips` rejects complemented matches but still computes the complemented reading.** If only that reading matches, the result is `EXCLUDED` with `exclusion_reason = 'strand_flip'`. This is a row in [ADR 0010](0010-multiallelic-matched-per-alt-bgen-stays-biallelic.md)'s exclusion table, not `NO_MATCH`. The variant exists, but policy excludes it. The log records that distinction.

Ambiguity depends on the complemented reading. Without that calculation, `is_ambiguous` cannot become true. The apparently cautious flag would silently score every palindromic SNP by assuming the forward strand. This would be more permissive than `--keep_ambiguous`, without an excluded row to expose the assumption. **Omitting the complemented calculation was therefore rejected**, despite its lower cost.

The flags remain independent. The matcher determines ambiguity before selecting a winning match. Without `--keep_ambiguous`, the matcher excludes ambiguous variants as `ambiguous`, regardless of the strand flag. Together, `--ignore_strand_flips --keep_ambiguous` retain an ambiguous variant initially, then exclude it if the complemented reading wins.

**The matching UDF makes the decision.** v2 computes all four flip strategies before applying `_label_flips`. That function sets `exclude=True` for any `match_type` containing `_flip` (`_match/label.py:296-310`). This separate step reflects the lazy polars pipeline, which cannot decide per variant.

The UDF already holds the whole target variant and returns one result, so it can decide directly. Collisions remain a separate SQL step after matching. A single UDF call cannot observe a collision between score variants.

**`--ignore_strand_flips` retains v2's spelling and default.** Flip matching remains enabled unless the user disables it. Its name describes the opposite action to the adjacent `--keep_*` flags. This is deliberate because existing v2 pipelines already use that name.

## Consequences

- **v3 differs from released v2 for multi-base alleles.** v2 uses a plain complement without a length guard. A pairing labelled `REFALT_FLIP` by v2 may be unmatched here, and a pairing unmatched by v2 may match through a flip here. Log comparisons will expose this difference. Single-base allele behaviour stays unchanged, so practical exposure is small.
- **The matcher leaves an indel's strand unresolved instead of resolving it incorrectly.** Opposite-strand indels become `NO_MATCH`. Match rates decrease, but the matcher avoids incorrect matches. The two allele strings alone cannot describe the required opposite-strand record.
- **`--ignore_strand_flips` can reduce match rates enough to reject a score.** `match_rate` counts only `match_summary = 'matched'`. Each excluded flip lowers the rate, which can fall below `min_overlap`. v2 behaves the same way because `skip_flip` also excludes the variant.
- **`match_flipped` counts zero for a variant excluded because of a flip.** The field uses `match_type ILIKE '%flip'`, but the type is now `EXCLUDED`. The summary's `strand_flip` count replaces that information exactly, without another column.
- **Nothing checks that all chromosome MATCH tasks receive identical flags.** Mixing `--ignore_strand_flips` settings within one sampleset produces an inconsistent log without an error. This limitation already exists for `--keep_ambiguous` and `--keep_multiallelic`. Nextflow supplies one parameter set to every task, so the limitation remains accepted.
