# 0009 — Only enforce the imputation sample floor when imputation is needed

**Refs:** ADR 0002 (consequences), ADR 0003

## Why

`fill_missing_dosages` (`_dosage.py:323`) calls `calculate_mean_dosage` unconditionally, and `calculate_mean_dosage` raises `ValueError` when `n_samples < n_minimum_impute` (default 50) *before* checking whether anything is missing. A shard with no missing calls and 40 samples crashes for no reason.

Harmless today because `n_samples` is always the whole cohort. Still harmless under alpha2, since rectangular batches impute within ~100k samples. Worth fixing anyway: it is a one-line guard, and it is the kind of latent bug that becomes a support burden the first time someone runs a small test cohort.

## Scope

- Skip imputation entirely when the dosage array has no missing values; only then can the floor apply.
- Independent of everything else — grab any time.

## Acceptance

- A 10-sample shard with no missing calls scores without raising.
- A 10-sample shard *with* missing calls still raises, naming the floor.
