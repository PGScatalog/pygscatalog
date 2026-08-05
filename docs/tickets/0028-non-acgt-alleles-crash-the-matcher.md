# 0028 — A spanning deletion or symbolic allele crashes MATCH

**Refs:** ADR 0010. Found while designing multiallelic support, but present on `main`
independently of it.

## Why

`is_biallelic_variant_ambiguous` raises on any allele outside `ACGT`
(`lib/score/_matchvariants.py:371-375`). The docstring around the matcher treats this as
unreachable, on the belief that the bcftools filter admits only plain SNPs and indels.
It does not:

```
--- bcftools view -V other ---            (excludes nothing)
1  300  star       A  G,*
1  400  symbolic   A  G,<NON_REF>
1  500  star_only  A  *
--- bcftools view -V other --min-alleles 2 --max-alleles 2 ---   (the current filter)
1  100  biallelic  A  G
1  500  star_only  A  *                    <- survives
```

`--exclude-types other` filters neither `*` nor `<NON_REF>`. Only `--max-alleles 2` is
doing any work, and a biallelic `A/*` record passes it. So:

```
match_variant(..., ref='A', alts=['*'])          -> ValueError: *
match_variant(..., ref='A', alts=['<NON_REF>'])  -> ValueError: <NON_REF>
```

The raise happens inside a DuckDB UDF, so the failure surfaces as a query error rather
than as anything a user can act on. A GATK-called VCF with a spanning deletion at a
scoring-file position is enough to trigger it. The comment at `lib/cache/_vcf.py:237-238`
claiming `--exclude-types` "excludes symbolic alleles" is wrong and should go with it.

## Scope

Make an allele that no score can match yield `NO_MATCH` instead of raising: no ambiguity
check, no comparison, next allele.

Keep the `ValueError` for genuinely corrupt allele strings. `*` and `<NON_REF>` are legal
VCF, with defined meanings — rejecting them is not strictness, it is a crash on
conformant input.

Fix the `_vcf.py` comment to say what the filter does.

## Out of scope

Removing `--max-alleles 2` (ticket 0029). This ticket is the crash alone, and is worth
landing first — without it, 0029 turns a rare crash into a common one.

## Acceptance

- `match_variant` against `alts=['*']` and `alts=['<NON_REF>']` returns `NO_MATCH`
  rather than raising.
- A VCF holding an `A/*` record at a scoring-file position scores without error, and
  the record is reported unmatched.
- One test covers each of `*` and a symbolic allele, driving `match_variant` directly.
