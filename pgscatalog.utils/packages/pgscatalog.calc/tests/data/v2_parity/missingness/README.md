# Missing hard calls

The fixture has two chromosomes and 52 diploid samples.
Each chromosome contains two A/C SNPs.
Called donors are heterozygous at position 1000 and homozygous REF at position 2000.
Each variant has 50 called donors.
Imputed dosages are integers: 0, 1 and 2.
The fixture contains no fractional dosage inputs or examples.

The missing calls follow this pattern:

| Samples | Missing calls |
| --- | --- |
| S00 | All calls |
| S01 | Both chromosome 1 calls |
| S02 | Both chromosome 2 calls |
| All remaining samples | None |

PGS990010 uses ALT effects on chromosome 1 with weights 2 and 0.25.
It uses REF effects on chromosome 2 with weights -0.5 and 0.5.
Mean imputation gives every sample a score sum of 2.5 and dosage sum of 4.

v3 divides the score sum by the observed allele count:

| Samples | Observed alleles | v3 average |
| --- | --- | --- |
| S00 | 0 | Null |
| S01 | 4 | 0.625 |
| S02 | 4 | 0.625 |
| All remaining samples | 8 | 0.3125 |

PLINK's denominator includes imputed alleles, so its average is 0.3125 for every sample.
PLINK reports only observed copies in `NAMED_ALLELE_DOSAGE_SUM`, even when it imputes the weighted sum.
v3's dosage sum includes imputed copies.
These are differences in scoring policy, not output format.

The unchanged v2 pipeline enables mean imputation for this cohort size.
A separate diagnostic reference disables mean imputation.
Its sums are 0 for S00, 0.5 for S01, 2 for S02 and 2.5 for all remaining samples.
The tests retain both sets of raw results and check each policy independently.

## Rejected inputs

`all_uncalled.chr1.vcf.gz` has no called donors.
LOAD excludes this uncalled VCF site.
No variants then match, so the score fails the matching threshold.

`too_small.chr1.vcf.gz` has four samples, including missing calls.
It must fail v3's minimum sample count of 50 for imputation.
These files test rejection. They are not additional chromosome shards for the reference workflow.

## Regenerate inputs

Run `generate.py <output-directory>` with the project Python environment.
The generator never writes reference outputs.
