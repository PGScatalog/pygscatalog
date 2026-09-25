# Arithmetic and partitioning fixture

Four fully called diploid samples, six A/C variants on each of chromosomes 1 and
2, and four scores exercise accession routing and decimal arithmetic:

- PGS990001: opposing large decimal weights, small weights and cancellation.
- PGS990002: overlapping variants, different weights, and the opposite effect
  allele at the first variant.
- PGS990003: variants disjoint from both preceding scores.
- PGS990004: additive, dominant and recessive contributions with dosages 0, 1, 2.

The tests calculate independent Decimal expectations directly from committed GT
calls and weights. Batch size 2 crosses actual scoring chunks without thousands
of artificial variants. Partition tests reorder sample columns independently per
chromosome, reverse chromosome/score input order, and split complete sample
batches. They reuse the unchanged reference outputs.

## Agreement and intentional differences

Selected alleles and weighted sums agree with the saved v2.3.0 run within the
precision of its original text outputs. PLINK writes six significant figures;
the tests sum half-unit rounding bounds for each original shard **before**
cancellation. v3's unrounded chromosome sums also face a much tighter independent
Decimal check, and its final output allows half a unit at six decimal places.

v2's multi-score files report shared variant-union allele counts and dosage sums,
including variants with zero weight for a particular accession. Those statistics
are not per-score statistics. Its final aggregation also includes denominators
from shards that do not contain that accession: every score has DENOM 36 here.
The corresponding v3 per-score observed-allele denominators are 16, 12, 4 and 12.
The tests explicitly preserve this difference instead of asserting average-score
parity or reproducing v2's denominator behaviour. Raw PGS990004 statistics agree:
each of its effect-type shards contributes one scored variant.

Normal tests launch no external scoring process. Run `generate.py` with this
directory as its argument to reproduce inputs; it never modifies references.
Fractional dosage inputs, missingness and sex chromosomes are absent here.
