# Zero weights preserve score membership

Three diploid samples have ALT dosages of 0/2, 1/1 and 2/0 across two A/C SNPs.
PGS990040 gives these loci weights of 0 and 2.
PGS990041 gives both loci a weight of zero.

Both scores match both loci.
Every sample therefore has four observed alleles and a dosage sum of two.
PGS990040 sums are 4, 2 and 0.
PGS990041 sums and averages are all zero.

An explicit zero coefficient still identifies a locus that belongs to the score.
It must remain distinct from a locus absent from that score when the code reshapes the weight matrix.
The test checks this membership policy from ADR 0015 and agreement with the committed raw v2 reference.

Run `generate.py <output-directory>` with the project Python environment to reproduce inputs.
The generator never writes or replaces reference files.
