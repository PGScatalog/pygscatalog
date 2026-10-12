# Imputation pools allele frequency across ploidy; the imputed dosage is per sample

## Status

Accepted. Follows from [ADR 0007](0007-ploidy-is-observed-per-call.md).

## Context

Missing dosages are imputed from the effect allele frequency observed in the shard. Both halves of that calculation assume every call is diploid:

```python
n_called_genotypes = da.sum(~da.isnan(dosage_array), axis=1)
n_called_alleles = n_called_genotypes * ZARR_PLOIDY        # _impute.py:51
effect_allele_freq = n_effect_alleles / n_called_alleles
mean_effect_allele_dosage = effect_allele_freq * ZARR_PLOIDY  # _impute.py:61
```

The result is a per-variant scalar, broadcast across every sample at `mean_effect_dosage[:, np.newaxis]` (`lib/score/_dosage.py:352`). `_impute.py:51` already carries the note `# replace this with a ploidy mask when X / Y / MT support added`.

Once ploidy is observed per call, both uses of `ZARR_PLOIDY` are wrong on X, Y and MT, and in different ways. The first inflates the denominator: a shard of 500,000 X calls from men contributes 500,000 alleles, not a million, so the frequency is halved. The second is wrong per sample rather than per variant — a man's missing X call should be imputed as one expected copy, a woman's as two.

That second correction changes the *shape* of the calculation, not just a constant. The imputed value stops being a per-variant scalar and becomes a per-`(variant, sample)` quantity.

A second question is whether the frequency itself should be pooled across ploidy or estimated separately for hemizygous and diploid samples. This can be asked without a sex column, because ploidy is observed: the strata are simply `ploidy == 1` and `ploidy == 2` at that variant.

## Decision

The frequency is **pooled** across ploidy — total effect alleles over total called alleles, which is a coherent allele frequency whatever mix of ploidies contributed to it. The **imputed dosage is per sample**: that sample's observed ploidy multiplied by the pooled frequency.

```python
ploidy = (gt != NO_ALLELE).sum(axis=2)                          # (n_variants, n_samples)
n_called_alleles = da.sum(da.where(~is_nan, ploidy, 0), axis=1)  # (n_variants,)
effect_allele_freq = n_effect_alleles / n_called_alleles
imputed = ploidy * effect_allele_freq[:, None]                   # (n_variants, n_samples)
filled = da.where(is_nan, imputed, dosage_array)
```

Nullploid samples drop out of both sides for free: they contribute 0 to the denominator and receive an imputed dosage of `0 × freq`, so a woman is never assigned a chrY genotype.

**Stratifying by ploidy was rejected.** It is robust to genuinely sex-differential X allele frequency, and it is free of any sex declaration. But it halves the sample size behind each estimate for a second-order effect, next to the much larger approximation of imputing at all; it is not what plink does, so it would make reproducing a published chrX score harder rather than easier; and it interacts badly with the `n_minimum_impute` floor, which would have to be satisfied within each stratum — a cohort of 1,000 with 12 men has a hemizygous stratum far below the default of 50, with no good answer for what to do then.

**Refusing to impute mixed-ploidy variants at all** was also rejected. It is the loudest option and fabricates nothing, but it discards real data and would leave X match rates worse than autosomal ones for a reason no user could act on.

## Consequences

- `fill_missing_dosages` (`lib/score/_dosage.py:305`) needs the ploidy array, so it can no longer take the dosage array alone. `calculate_mean_dosage` likewise. This is the interface change; the arithmetic is small.
- The `n_minimum_impute` floor keeps its current meaning — a count of samples, checked once — because pooling means there is only one estimate. Had stratification been chosen, the floor would have needed redefining as a count of alleles or of samples per stratum.
- Pooling assumes the effect allele frequency does not differ between hemizygous and diploid samples at a variant. Where that assumption fails — sex-differential selection, or a cohort whose men and women differ in ancestry composition — imputed X dosages are biased. Only imputed values are affected; observed calls are untouched. The exposure is therefore bounded by the X missingness rate, which is the quantity to look at before revisiting this.
- Score values change for any accession carrying X variants that was previously being scored with those variants dropped. See [ADR 0007](0007-ploidy-is-observed-per-call.md) for why the `min_overlap` gate itself cannot move against a user.
- A NaN frequency is still possible where a variant has no called alleles at all, and is still caught by the existing post-write NaN check in `write_dosage_arrays` (`lib/score/_dosage.py:226`). Nullploidy makes this reachable in a new way — an all-female shard at a chrY variant has zero called alleles by construction — so the check's error message should distinguish "nothing called here" from "nobody can carry this", or the first user to hit it will read it as data corruption. **Implemented:** `no_called_alleles_reason` discriminates the two on `all(ploidy == 0, axis=1)` over the variants that went NaN.
- **The zero-denominator guard had to stop evaluating the division.** The original `da.where(n_called > 0, n_effect / n_called, nan)` computes a real `0/0` in the rejected branch and emits `RuntimeWarning: invalid value encountered in divide`, which this suite escalates to an error. It is now `n_effect / da.where(n_called > 0, n_called, nan)` — same NaN result, no warning. Latent before this ADR only because nothing exercised a zero denominator; nullploidy makes it ordinary.
- **This decision is currently undermined on the VCF path, and not by its own arithmetic.** The imputed dosage is `ploidy × freq`, which is only as good as the ploidy recorded for a *missing* call — and `vcf_get_genotypes` stores every wholly missing call as diploid whatever its `GT` length, departing from [ADR 0007](0007-ploidy-is-observed-per-call.md)'s `hemizygous no-call [255, 254]`. So a man's missing chrX call imputes `2 × freq` while contributing 0 alleles to his `allele_count`: his score is inflated and nothing warns. The pooled *frequency* is correct — it is computed from observed calls, whose ploidy is read correctly — so only imputed values for hemizygous no-calls are affected. bgen is unaffected, because `BgenVar.ploidy` stays correct for missing calls. The per-sample half of this decision is therefore fully exercised only through bgen until that encoding is fixed.
