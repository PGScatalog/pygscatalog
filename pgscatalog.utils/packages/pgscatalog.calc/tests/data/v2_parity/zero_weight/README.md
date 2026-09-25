# Zero weights preserve score membership

Three diploid samples have opposite ALT dosages at two A/C SNPs: 0/2, 1/1 and
2/0. PGS990040 weights these loci by 0 and 2; PGS990041 weights both by zero.
Both scores match both loci, so every sample has four observed alleles and a
dosage sum of two. PGS990040 sums are 4, 2 and 0; PGS990041 sums and averages
are all zero. An explicit zero coefficient must not be confused with a locus
that is absent from a score after reshaping the weight matrix.

This tests the matched-row membership policy in ADR 0015 and agreement with
the committed raw v2 reference. Run `generate.py <output-directory>` using
the project Python environment to regenerate inputs; reference files are
never generated or overwritten by that script.
