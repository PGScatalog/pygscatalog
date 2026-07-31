# Chromosome names are canonicalised at the boundary; `XY` collapses into `X`

## Status

Accepted. Enables [ADR 0007](0007-ploidy-is-observed-per-call.md).

## Context

Chromosome names are load-bearing in two places that must agree: the zarr hierarchy is `/{sampleset}/{chromosome}/{filename}/genotypes` ([ADR 0003](0003-shards-are-zarr-leaf-groups.md)), and matching joins scoring file variants to target variants on `(chr_name, chr_pos)`. While everything was 1–22 the dialects happened to line up. X, Y and MT break that.

The three vocabularies in play are all different:

- **PGS Catalog scoring files** use `X`, `XY`, `MT`. This is not inferred — `add_chrom_padding`'s docstring names those exact values as the non-numeric cases it must pass through (`lib/cache/_bgen.py:404`).
- **Target genomes** use `X` or `chrX`, `MT` or `chrM` or `M`, with or without zero padding.
- **plink**, a format the reader interface is expected to grow, uses numeric codes: 23 = X, 24 = Y, 25 = XY (pseudoautosomal), 26 = MT.

Handling today is asymmetric and only half exists. The bgen path normalises target names — `normalise_chrom` strips a `chr` prefix and leading zeros (`_bgen.py:183`), with `add_chrom_padding` and `add_chrom_prefix` translating back into the target's dialect. The VCF path normalises **nothing**: it passes scoring-file chromosome names straight to `bcftools view -R`, so a VCF using `chrX` returns no regions and reports the variants as unmatched rather than failing.

`XY` is the sharp edge, and it is purely a naming problem. Under ADR 0007 a pseudoautosomal variant needs no biological special-casing — it is diploid in the file because it is diploid, and ploidy is observed. But a scoring file variant labelled `XY` still has to *match* a target genome that calls that locus `X`, and no amount of correct ploidy handling makes two different strings join.

## Decision

There is one internal chromosome vocabulary: **`1`–`22`, `X`, `Y`, `MT`**. Names are canonicalised **once, at the boundary** — every `GenomeFileHandler` on the way in, and the scoring file loader on the way in — and the Cache stores canonical names only.

```
scoring file   XY, X, MT      ->  X, X, MT
VCF            chrX, chrM     ->  X, MT
bgen           0X, chrX       ->  X
plink          23, 25, 26     ->  X, X, MT
```

`XY` collapses into `X`. `M` and `chrM` collapse into `MT`. `VALID_CHROMOSOMES` becomes `{1..22, X, Y, MT}` and stops being a list of what is implemented.

This is the same principle [ADR 0003](0003-shards-are-zarr-leaf-groups.md) already applies to genotypes: after LOAD, no stage knows what the input looked like. Extending it to chromosome names means MATCH never compares across dialects, and the zarr hierarchy stays legible to a domain expert inspecting the store — `/{sampleset}/X/...` regardless of whether the input said `chrX`, `0X` or `23`.

**Keeping `XY` as a distinct internal chromosome was rejected.** It preserves a distinction plink models, but every target genome that calls the pseudoautosomal region `X` would then need its PAR variants relabelled, which requires exactly the build-specific PAR coordinate map ADR 0007 declined to ship. The distinction buys nothing once ploidy is observed rather than derived, because PAR's only significance was telling you where a man is diploid — and the file already says.

**Per-reader normalisation only** — extending the bgen approach and leaving the cache to hold whatever the target used — was rejected because it pushes dialect reconciliation into MATCH, which then has to compare `chrX` against `X` against `23` for every sampleset.

## Consequences

- The VCF path's total absence of normalisation is fixed as a side effect rather than as separate work. A `chr`-prefixed VCF silently matching nothing is a pre-existing bug for chromosomes 1–22 too, not something X introduced.
- `vcf_buffer_variants` sorts regions with `key=lambda t: (int(t[0]), t[1])` (`lib/cache/_vcf.py:65`), which raises on any non-numeric canonical name. The sort key has to become dialect-aware; canonicalising does not save it, because `X` stays `X`.
- Canonical names are baked into the zarr hierarchy, so a cache written before this decision and one written after are not interchangeable. As with [ADR 0003](0003-shards-are-zarr-leaf-groups.md)'s depth change, alpha caches are discarded rather than migrated.
- A scoring file carrying both an `XY` variant and an `X` variant at the same position now produces two rows at one canonical position. **Measured rather than left open:** it is neither collapsed as a duplicate nor rejected as an error. `score_variant_table` holds both rows with their two different weights; `is_duplicated` stays `False` on both, because pgscatalog-format computes it before `XY` becomes `X` and nothing recomputes it afterwards. `allele_match_table` matches both against the same `target_row_nr` — the PK is `(sampleset, accession, row_nr, zarr_group)`, so nothing dedups them. `create_wide_weights_table`'s PIVOT then groups on a key identical for both rows and aggregates with `coalesce(any_value(effect_weight), 0)`, yielding **one weight-matrix row carrying one of the two weights, with which one unspecified** — both tie on the insert's `ORDER BY`. So the `GROUP BY` prevents a double count in the score, but a weight is discarded silently and the `min_overlap` numerator counts one locus twice. `is_duplicated` is therefore the wrong place to look; a `(accession, chr_name, chr_position)` uniqueness check on `score_variant_table` is the natural fix if the answer is "reject". Any two variants that canonicalise to one position hit this — `XY`/`X` is just the first way to reach it. Worth noting that `_weight_matrix.py` states the broken assumption in a comment beside the aggregation — *"the pivoted values (effect weights) are unique on accession, any_value is a workaround, it's no-op"* — which was true before canonicalisation and is what this decision invalidates.
- Collapsing `XY` into `X` is not reversible from the cache. If a future decision needs PAR as a distinct region — a PAR-aware quality check, say — it must come from a coordinate map, not from a chromosome name.
