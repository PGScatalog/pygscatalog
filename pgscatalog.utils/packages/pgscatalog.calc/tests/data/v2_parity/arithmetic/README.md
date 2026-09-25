# Arithmetic and partitioning fixture

The fixture has four diploid samples with complete genotype calls.
Chromosomes 1 and 2 each contain six A/C variants.
Four scores test the assignment of weights to accessions and decimal arithmetic:

- PGS990001 uses opposing large decimal weights, small weights and cancellation.
- PGS990002 uses overlapping variants with different weights. It uses the opposite effect allele at the first variant.
- PGS990003 uses variants that occur in neither of the preceding scores.
- PGS990004 uses additive, dominant and recessive contributions with dosages 0, 1 and 2.

The tests calculate independent expectations from committed GT calls and weights with `Decimal`.
A batch size of 2 tests calculations across actual scoring chunks without thousands of artificial variants.
Partition tests change the sample-column order separately for each chromosome.
They reverse the chromosome and scoring-file input order.
They also divide the samples into complete batches.
All these tests use the same reference outputs.

## Agreement and intentional differences

Selected alleles and weighted sums agree with the saved v2.3.0 run within the precision of its original text outputs.
PLINK writes six significant figures.
The tests add half-unit rounding bounds for each original shard before contributions cancel.
v3's unrounded chromosome sums also receive an independent `Decimal` check with a smaller error bound.
Checks of the final output allow half a unit at six decimal places.

v2 groups multiple scores into files with shared allele counts and dosage sums.
These statistics describe the union of the variants, including variants with zero weight for a particular accession.
They do not describe each score separately.
v2's final aggregate also includes denominators from shards without that accession.
Every score therefore has DENOM 36 in this reference.

v3 uses a separate observed-allele denominator for each score.
The denominators are 16, 12, 4 and 12, in accession order.
The tests retain this difference and do not require equal average scores.
The raw PGS990004 statistics agree because each file for an effect type contributes one scored variant.

Normal tests launch no external scoring process.
Run `generate.py` with this directory as its argument to reproduce inputs.
The generator never changes references.
This fixture excludes fractional dosage inputs, missingness and sex chromosomes.
