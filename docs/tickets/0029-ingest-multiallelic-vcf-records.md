# 0029 — Multiallelic VCF records never reach the cache

**Refs:** ADR 0010.

## Why

`--keep_multiallelic` is unreachable code. VCF ingest passes `--min-alleles 2
--max-alleles 2` to bcftools (`lib/cache/_vcf.py:262-265`), so a multiallelic record is
dropped before the cache exists and the flag has nothing to act on.

Worse than dropped: the position is then *misreported*. It never enters
`seen_positions`, so `add_missing_positions_to_lists` pads it with `ref=None`
(`lib/cache/targetvariants.py:303-310`), `update_match_table`'s
`WHERE target.ref IS NOT NULL` filters it out, and the match log calls it unmatched. The
user is told a variant was not found in their data. It was found, and deliberately
discarded.

Allele index headroom is not a constraint here — genotypes are `uint8` with 254 and 255
reserved (`lib/constants.py:9-18`), so a record can carry 253 alternative alleles before
the encoding breaks. Nothing in the dosage path needs changing either: dosage counts
calls equal to `effect_allele_idx` (`lib/score/_dosage.py:317`), which is already
index-based and allele-count agnostic.

## Scope

Drop `--max-alleles 2` so multiallelic records enter the cache with their full allele
list.

Alleles are positional: a genotype call names the allele it carries by index, so an
allele can never be removed from a record without renaming every allele after it.
Dropping `*` from `[C,*,G]` renumbers `G` from 3 to 2 and every dosage for it becomes the
wrong allele. `TargetVariants.alts` therefore holds every alternative allele the record
declared, in the order it declared them.

## Out of scope

Matching against the second and later alternative alleles (ticket 0030). After this
ticket the records are cached and `--keep_multiallelic` becomes reachable, but the
matcher still examines only `alts[0]`.

bgen (ticket 0032).

## Acceptance

- A VCF holding `A -> C,G` at a scoring-file position caches one variant with
  `ref='A'`, `alts=['C','G']`, and genotype indices 0, 1 and 2 preserved as written.
- That position no longer appears in the padded-absent set, and no longer reports as
  unmatched when it was assessed.
- A round-trip test reads a multiallelic record out of the cache and confirms a
  `1/2` call is stored as `[1, 2]`.
