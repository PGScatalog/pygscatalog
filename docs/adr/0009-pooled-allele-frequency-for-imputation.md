# Imputation pools allele frequency across ploidy and calculates dosage per sample

## Status

Accepted. Follows from [ADR 0007](0007-ploidy-is-observed-per-call.md).
[ADR 0016](0016-pool-imputation-frequencies-across-the-sampleset.md) extends frequency estimation from one shard to the whole sampleset; pooling across ploidy remains unchanged.

## Context

Imputation uses the effect allele frequency observed in the shard. At this decision, both parts of the calculation assume diploid calls:

```python
n_called_genotypes = da.sum(~da.isnan(dosage_array), axis=1)
n_called_alleles = n_called_genotypes * ZARR_PLOIDY        # _impute.py:51
effect_allele_freq = n_effect_alleles / n_called_alleles
mean_effect_allele_dosage = effect_allele_freq * ZARR_PLOIDY  # _impute.py:61
```

The result is one scalar per variant. `mean_effect_dosage[:, np.newaxis]` broadcasts it across samples (`lib/score/_dosage.py:352`). `_impute.py:51` already contains the note `# replace this with a ploidy mask when X / Y / MT support added`.

Both uses of `ZARR_PLOIDY` become incorrect on X, Y and MT when readers observe ploidy per call. The first inflates the denominator. A shard with 500,000 X calls from men contributes 500,000 alleles, not a million. The incorrect denominator halves the frequency.

The second error affects each sample. A man's missing X call needs one expected copy, while a woman's needs two. This correction changes the calculation's shape, not only a constant. The imputed value becomes a per-`(variant, sample)` quantity instead of a scalar per variant.

The frequency could pool ploidies or use separate estimates for hemizygous and diploid samples. Observed ploidy permits either choice without a sex column. The strata would be `ploidy == 1` and `ploidy == 2` at each variant.

## Decision

**Pool the frequency across ploidies.** Divide total effect alleles by total called alleles. This gives a coherent allele frequency for any mixture of ploidies. **Calculate imputed dosage per sample** as its observed ploidy multiplied by the pooled frequency.

```python
ploidy = (gt != NO_ALLELE).sum(axis=2)                          # (n_variants, n_samples)
n_called_alleles = da.sum(da.where(~is_nan, ploidy, 0), axis=1)  # (n_variants,)
effect_allele_freq = n_effect_alleles / n_called_alleles
imputed = ploidy * effect_allele_freq[:, None]                   # (n_variants, n_samples)
filled = da.where(is_nan, imputed, dosage_array)
```

Nullploid samples contribute 0 to the denominator and receive an imputed dosage of `0 × freq`. Thus, imputation never assigns a chrY genotype to a woman.

**This ADR rejects stratification by ploidy.** It handles real sex differences in X allele frequency without a sex declaration. However, it halves the samples behind each estimate to address a secondary effect. Imputation itself introduces a much larger approximation.

Stratification also differs from plink and would complicate reproduction of published chrX scores. Each stratum would need to satisfy `n_minimum_impute`. A cohort of 1,000 with 12 men has a hemizygous stratum far below the default of 50. There is no satisfactory policy for that case.

**This ADR also rejects refusing to impute mixed-ploidy variants.** It would report the limitation clearly and fabricate nothing. However, it would discard real data and reduce X match rates relative to autosomes. Users could not correct the cause.

## Consequences

- `fill_missing_dosages` needs the ploidy array in addition to the dosage array (`lib/score/_dosage.py:305`). `calculate_mean_dosage` needs the same change. The interface changes, but the arithmetic remains small.
- `n_minimum_impute` remains a sample count, checked once, because pooling gives one estimate. Stratification would require an allele count or a sample count for each stratum instead.
- Pooling assumes equal effect allele frequencies between hemizygous and diploid samples at a variant. Sex-specific selection or ancestry differences between men and women can violate this assumption and bias imputed X dosages. Observed calls remain unchanged. The X missingness rate limits exposure and should inform any reconsideration.
- Scores change for accessions with X variants that the previous calculation discarded. [ADR 0007](0007-ploidy-is-observed-per-call.md) explains why the `min_overlap` gate cannot become less permissive.

### No called alleles

A variant with no called alleles can still produce a NaN frequency. The existing check after writing detects it in `write_dosage_arrays` (`lib/score/_dosage.py:226`). Nullploidy introduces another route: an all-female shard has zero called alleles at a chrY variant.

The error should distinguish "nothing called here" from "nobody can carry this". Otherwise, users could mistake the latter for data corruption. **Implemented:** `no_called_alleles_reason` checks `all(ploidy == 0, axis=1)` over variants with NaN values.

**The zero-denominator guard needed to avoid division itself.** The original `da.where(n_called > 0, n_effect / n_called, nan)` evaluates `0/0` in the rejected branch. It emits `RuntimeWarning: invalid value encountered in divide`, which the test suite treats as an error.

The expression is now `n_effect / da.where(n_called > 0, n_called, nan)`. It produces the same NaN result without a warning. The defect previously remained hidden because no test exercised a zero denominator. Nullploidy makes that case ordinary.

### Recorded VCF limitation

This section records the VCF limitation reported at the time. It concerns the stored ploidy of missing calls, not the imputation arithmetic. The imputed dosage is `ploidy × freq`. At that point, `vcf_get_genotypes` stored every wholly missing call as diploid, regardless of `GT` length. This contradicted [ADR 0007](0007-ploidy-is-observed-per-call.md)'s `hemizygous no-call [255, 254]` encoding.

A man's missing chrX call therefore imputed `2 × freq` while contributing 0 to `allele_count`. The score increased incorrectly without a warning. The pooled frequency remained correct because observed calls retained their correct ploidy. The defect affected only imputed hemizygous no-calls.

The defect did not affect bgen because `BgenVar.ploidy` remained correct for missing calls. Until the VCF encoding changed, only bgen fully exercised the per-sample part of this decision.
