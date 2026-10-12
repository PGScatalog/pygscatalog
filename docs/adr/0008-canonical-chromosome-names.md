# Chromosome names become canonical at the boundary. `XY` becomes `X`

## Status

Accepted. Enables [ADR 0007](0007-ploidy-is-observed-per-call.md).

## Context

Chromosome names must agree between the zarr hierarchy and matching. The hierarchy is `/{sampleset}/{chromosome}/{filename}/genotypes` ([ADR 0003](0003-shards-are-zarr-leaf-groups.md)). Matching joins scoring variants to target variants on `(chr_name, chr_pos)`. The naming conventions agreed while every chromosome was 1–22. X, Y and MT break that agreement.

The three sources use different names:

- **PGS Catalog scoring files** use `X`, `XY`, `MT`. The `add_chrom_padding` docstring explicitly lists these non-numeric values (`lib/cache/_bgen.py:404`). This is evidence, not an inference.
- **Target genomes** use `X` or `chrX`, and `MT`, `chrM` or `M`. They can also use zero padding.
- **plink** uses numeric codes: 23 = X, 24 = Y, 25 = XY (pseudoautosomal), 26 = MT. The reader interface is expected to support plink later.

At this decision, only the bgen path normalises target names. `normalise_chrom` removes a `chr` prefix and leading zeros (`_bgen.py:183`). `add_chrom_padding` and `add_chrom_prefix` restore the target's naming convention.

The VCF path normalises nothing. It passes scoring-file chromosome names directly to `bcftools view -R`. A VCF using `chrX` returns no regions and reports unmatched variants instead of an error.

`XY` presents a naming problem. Under ADR 0007, pseudoautosomal variants need no separate biological handling because the file records their diploidy. A score variant labelled `XY` must still match a target locus labelled `X`. Correct ploidy handling cannot make different strings join.

## Decision

The internal chromosome names are **`1`–`22`, `X`, `Y`, `MT`**. Every `GenomeFileHandler` and the scoring file loader canonicalise names once, when they read inputs. The Cache stores only canonical names.

```
scoring file   XY, X, MT      ->  X, X, MT
VCF            chrX, chrM     ->  X, MT
bgen           0X, chrX       ->  X
plink          23, 25, 26     ->  X, X, MT
```

`XY` becomes `X`. `M` and `chrM` become `MT`. `VALID_CHROMOSOMES` becomes `{1..22, X, Y, MT}` rather than a list of implemented chromosomes.

[ADR 0003](0003-shards-are-zarr-leaf-groups.md) already applies this principle to genotypes. After LOAD, no stage needs to know the input representation. MATCH therefore compares only canonical names. A domain expert sees `/{sampleset}/X/...` in the store whether the input used `chrX`, `0X` or `23`.

**This ADR rejects a distinct internal `XY` chromosome.** This preserves a distinction in plink, but requires relabelling PAR variants in targets that use `X`. Relabelling needs the build-specific PAR coordinate map that ADR 0007 rejected. The distinction adds no information when the file already records where a man is diploid.

**This ADR rejects normalisation only within readers.** This would extend the bgen approach and let the Cache retain target names. MATCH would then need to reconcile `chrX`, `X` and `23` for every sampleset.

## Consequences

- This also corrects the VCF path's missing normalisation. A `chr`-prefixed VCF already silently matched nothing on chromosomes 1–22. X did not introduce that defect.
- `vcf_buffer_variants` uses `key=lambda t: (int(t[0]), t[1])` to sort regions (`lib/cache/_vcf.py:65`). This fails for non-numeric canonical names. The sort key must handle those names because canonicalisation leaves `X` as `X`.
- The zarr hierarchy contains canonical names. Caches from before and after this decision are not interchangeable. As with [ADR 0003](0003-shards-are-zarr-leaf-groups.md)'s depth change, discard alpha caches instead of migrating them.

### Recorded collision behaviour

Measurements during this decision established the following behaviour. A scoring file with `XY` and `X` variants at the same position produces two rows at one canonical position. The implementation neither merges them as duplicates nor rejects them.

`score_variant_table` retains both rows and their different weights. Both retain `is_duplicated = False`. pgscatalog-format computes this flag before `XY` becomes `X`, and nothing recomputes it afterwards.

`allele_match_table` matches both rows to the same `target_row_nr`. Its primary key is `(sampleset, accession, row_nr, zarr_group)`, so it removes neither row. The PIVOT in `create_wide_weights_table` groups both rows under the same key. It applies `coalesce(any_value(effect_weight), 0)` and produces one weight-matrix row with one unspecified weight. Both rows tie under the insert's `ORDER BY`.

The `GROUP BY` prevents a double score contribution, but silently discards one weight. The `min_overlap` numerator counts the locus twice. Checking `is_duplicated` cannot detect this problem. If the decision is to reject it, the natural check is uniqueness of `(accession, chr_name, chr_position)` in `score_variant_table`.

Any two variants that canonicalise to one position can cause this problem. `XY` and `X` provide the first example. A comment beside the aggregation in `_weight_matrix.py` states the invalid assumption: *"the pivoted values (effect weights) are unique on accession, any_value is a workaround, it's no-op"*. That assumption held before canonicalisation, which invalidates it.

Collapsing `XY` into `X` cannot be reversed from the Cache. Any future PAR-specific decision, such as a quality check, must use a coordinate map rather than a chromosome name.
