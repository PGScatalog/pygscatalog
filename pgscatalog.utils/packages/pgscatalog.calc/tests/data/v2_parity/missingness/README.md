# Missing hard calls

The two chromosomes each contain two A/C SNPs and 52 diploid samples. Donors
are heterozygous at position 1000 and homozygous REF at position 2000. S00 has
all calls missing, S01 lacks both chromosome 1 calls and S02 lacks both
chromosome 2 calls. Each variant has 50 called donors. Integer imputed dosages
cover 0, 1 and 2; there are no fractional dosage inputs or examples.

PGS990010 uses ALT effects on chromosome 1 with weights 2 and 0.25, and REF
effects on chromosome 2 with weights -0.5 and 0.5. Mean imputation gives every
sample a sum of 2.5 and dosage sum of 4. Observed allele counts are 0, 4, 4 and
8 for S00, S01, S02 and all remaining samples. v3 divides by these observed
counts: the corresponding averages are null, 0.625, 0.625 and 0.3125. PLINK's
mean-imputed denominator includes missing alleles, so its average is 0.3125
for everyone. PLINK reports only observed copies in `NAMED_ALLELE_DOSAGE_SUM`,
even when it imputes the weighted sum; v3's dosage sum includes imputed copies.
These are explicit policy differences, not output formatting comparisons.

The unchanged v2 pipeline enables mean imputation for this cohort size. A
separate diagnostic reference explicitly disables mean imputation: its sums
are 0, 0.5, 2 and 2.5 respectively. The tests retain both sets of raw results
and verify each policy independently.

`all_uncalled.chr1.vcf.gz` has no called donors. LOAD excludes this uncalled
VCF site, so no variants match and the score fails the matching threshold.
`too_small.chr1.vcf.gz` has four samples, including missing calls, and must
fail v3's imputation floor of 50 samples. These are rejection inputs, not
additional chromosome shards for the reference workflow.

Regenerate inputs with `generate.py <output-directory>` using the project
Python environment. The generator never writes reference outputs.
