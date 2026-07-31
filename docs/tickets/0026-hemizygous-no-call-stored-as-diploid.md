# 0026 — A hemizygous no-call is stored as diploid missing, and inflates the score

**Refs:** ADR 0007. Found while implementing 0020 and 0021.

## Why

ADR 0007's encoding table lists five cases, and the VCF reader implements four of them:

```
hemizygous no-call   [255, 254]
```

`vcf_get_genotypes` (`lib/cache/_vcf.py`) collapses *any* wholly missing call to
`MISSING_GENOTYPE_TUPLE` — `[255, 255]` — whatever its `GT` length. So a man's missing
chrX call is stored as a missing **diploid** call.

The consequence is a wrong score rather than a wrong shape. 0021 makes imputation
ploidy-aware, so the imputed dosage is `ploidy × freq`; a call stored at ploidy 2
imputes `2 × freq` where the man can only carry one allele. His `allele_count`, derived
from the same array, counts the variant at ploidy 2 as well. The score is inflated and
nothing warns.

The reason recorded in the docstring is that "VCF writes `.` for both a haploid no-call
and a diploid one". That is not what pysam does — the two are distinguishable:

```
./.  ->  GT=(None, None)  len=2
.    ->  GT=(None,)       len=1
```

The claim came from this repo's own fixtures, which are built with
`sample["GT"] = None`; pysam emits `.` for that, so honouring `GT` length would have
made every fixture's autosomal missing call look hemizygous. A conformant VCF writes
`./.` for a diploid no-call, so on real input the length is the right signal.

## Scope

Decide whether `GT` length is trusted for wholly missing calls, and make the reader
agree with ADR 0007's table either way.

The reason this is a ticket rather than a one-line fix: trusting length means trusting
every upstream writer to be conformant, and a writer that emits `.` for a diploid
no-call would silently halve those samples' imputed dosage — the same class of error in
the other direction. Whether that is acceptable, or whether a mixed-ploidy variant with
missing calls should be rejected outright (0020's original position), is a decision this
ticket has to make rather than assume.

Fix the fixtures so they stop being the reason for the answer.

## Out of scope

bgen, which already gets this right — `BgenVar.ploidy` reports the true ploidy of a
sample whose call is missing, so `bgen_to_hardcall` writes `[255, 254]` correctly.

## Acceptance

- A man's missing chrX call in a conformant VCF imputes to `1 × freq`, not `2 × freq`,
  and counts 1 toward his `allele_count`. This is ticket 0021's first acceptance
  criterion, which is currently only reachable by driving `write_dosage_arrays`
  directly.
- Whatever is decided for a `.` written by a non-conformant diploid writer, it is
  decided explicitly and the decision is recorded — not left to fixture construction.
- The `vcf_get_genotypes` docstring states what the code does.
