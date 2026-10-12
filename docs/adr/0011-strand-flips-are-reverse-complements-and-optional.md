# A strand flip is a reverse complement, and matching one is optional

## Status

Accepted. Defines the complemented reading that
[ADR 0010](0010-multiallelic-matched-per-alt-bgen-stays-biallelic.md) computes
ambiguity from, and adds the flag that declines to match on it.

## Context

A scoring file and a target genome can report the same variant from opposite strands.
The strand-flipped reading of a sequence is its **reverse** complement: a sequence is
read 5'→3' on both strands, so flipping strand reverses the order of the bases as
well as complementing each one.

```
'A'    complement 'T'     reverse complement 'T'    <- agree
'AT'   complement 'TA'    reverse complement 'AT'   <- differ
'GC'   complement 'CG'    reverse complement 'GC'
'ATG'  complement 'TAC'   reverse complement 'CAT'
```

For a single base the reversal is invisible, which is why a plain complement passes
casual inspection: SNPs are nearly all of the exposure. For anything longer the two
disagree, and a plain complement matches sequences a strand flip does not produce
while missing the ones it does — a score variant `AC/CC` would be `REFALT_FLIP`
against a target `TG/GG`, which is not what the other strand of `AC/CC` carries, and
unmatched against `GT/GG`, which is.

Indels are the reason this needs a decision rather than a definition. A VCF indel
allele carries a shared anchor base at the start (`REF=A ALT=AT`), and that anchor is
*positional* — it stays leftmost on whichever strand the record is read.

Two consumers depend on the answer, not one. The complemented reading decides
`REFALT_FLIP` / `ALTREF_FLIP` matches, and it is also what
[ADR 0010](0010-multiallelic-matched-per-alt-bgen-stays-biallelic.md) computes
ambiguity from: a score variant is ambiguous when a complemented reading of it also
matches an allele the target carries.

plink2 contributes no strand logic to compare against. v2's flipped allele becomes
`matched_effect_allele` and is renamed `effect_allele` in the plink2 scorefile
(`_match/plink.py:96`), so plink2 matches a literal string against the pvar and has
no opinion about strand. There is no external behaviour to stay compatible with: the
decision is entirely ours.

Some users know which strand their target genotypes are on, and would rather lose a
variant than match one on an inferred strand. v2 gives them `--ignore_strand_flips`
(`match_cli.py:246`), so v3 needs an answer for them too.

## Decision

**A complemented reading is a reverse complement.** `reverse_complement` complements
each base and reverses the sequence. The SQL macro used to derive
`matched_effect_allele` for the match log is the same function, so the log never
names a different allele from the match that produced it.

**No complemented reading is attempted where the two alleles differ in length.** An
indel can still match forward (`REFALT`, `ALTREF`, `REF_NO_OA`, `ALT_NO_OA`); it
simply never matches via a complemented reading, and is therefore never ambiguous
either.

The rejected alternative is reverse-complementing indels anyway. Because the anchor
base is positional, the strand-flipped record sits at a different position with
different alleles, and no string operation on the two alleles in hand recovers it.
Reverse-complementing an indel produces something plausible-looking and wrong, which
is exactly the failure mode this ADR exists to remove. Declining is the honest
answer: an indel written on the other strand is a record this matcher cannot
recognise, and saying so is better than matching the wrong one.

**The length comparison is made between the target's reference allele and the one
alternative allele being read**, and it gates the whole complemented reading rather
than one branch of it. "The two alleles are the same length" is a property of the
*pairing* being read, never of one allele alone: `A` says nothing about whether the
site is an indel, and a score variant with no other allele has no second allele of
its own to compare with. The target pair is the one that is always available —
biallelic or multiallelic, with an other allele or without — and it is per
alternative allele, so a multiallelic variant carrying both a SNP and an indel
alternative gets the right answer for each.

It is also sufficient for the readings that *do* have a score-side pair, because a
reverse complement preserves length: a complemented match against an equal-length
target pair forces the score alleles to be equal-length too. Checking the score side
as well would add a condition that can never fire.

**`--ignore_strand_flips` declines to match on a complemented reading, and the
reading is still always computed.** A score variant whose only matching reading was a
complemented one is `EXCLUDED` with `exclusion_reason = 'strand_flip'` — a row in
[ADR 0010](0010-multiallelic-matched-per-alt-bgen-stays-biallelic.md)'s exclusion
table, not a `NO_MATCH`. The variant was found; policy is what discarded it, and the
log says so.

The reading is still computed because ambiguity is *defined* by it. Skip the
computation and `is_ambiguous` is structurally unreachable, so a flag whose name
promises caution would silently score every palindromic SNP on a forward-strand
assumption — strictly more permissive than `--keep_ambiguous`, with no excluded row
to warn anyone. **Not attempting the reading at all is therefore rejected**, cheaper
though it is.

The two flags stay orthogonal as a result. Ambiguity is decided first, before a
winning match is known, so an ambiguous variant is excluded as `ambiguous` whatever
the strand flag says; `--ignore_strand_flips --keep_ambiguous` together keep an
ambiguous variant and then exclude it if the reading that won was complemented.

**The decision is made inside the matching UDF, not as a post-match pass.** v2
computes all four flip strategies, then rescues the result afterwards with
`_label_flips` setting `exclude=True` on any `match_type` containing `_flip`
(`_match/label.py:296-310`). That shape is an artefact of a lazy polars pipeline that
cannot decide per variant. The UDF already holds the whole target variant and yields
one result, so it decides there. Collisions remain a post-match SQL pass, because a
collision is the one exclusion that is invisible inside a single UDF call.

**`--ignore_strand_flips` keeps v2's spelling and default.** Flip matching stays on
unless asked. The name reads backwards from the `--keep_*` flags beside it, and that
is deliberate: users arrive from v2 with the name in their pipelines.

## Consequences

- **v3 diverges from released v2 on multi-base alleles.** v2 uses a plain complement
  with no length guard, so a multi-base pairing it reports `REFALT_FLIP` may be
  unmatched here, and one it reports unmatched may be a flipped match here. Anyone
  diffing v3 match logs against v2 will see it. Single-base alleles are unaffected,
  so the practical exposure is small.
- **An indel's strand is unresolvable rather than wrongly resolved.** Scoring files
  carrying indels written on the opposite strand from the target lose those variants
  to `NO_MATCH`. That costs match rate and buys correctness; there is no
  representation of a strand-flipped indel that the two alleles alone can reach.
- **`--ignore_strand_flips` costs match rate, and can cost a score.** `match_rate`
  sums only `match_summary = 'matched'`, so every flip-excluded variant lowers it and
  an accession can fall under `min_overlap` and go unscored. This matches v2, where
  `skip_flip` also excludes rather than rescues.
- **`match_flipped` counts zero for a flip-excluded variant.** It is derived from
  `match_type ILIKE '%flip'` and the type is now `EXCLUDED`. The `strand_flip`
  exclusion count in the summary log replaces it exactly, so no column is added.
- **Nothing checks that every chromosome's MATCH task was given the same flags.** A
  sampleset matched with `--ignore_strand_flips` on some chromosomes and not others
  has an incoherent match log and no error. This is pre-existing with
  `--keep_ambiguous` and `--keep_multiallelic`; Nextflow passes one parameter set to
  every task, which is why it is tolerated rather than solved.
