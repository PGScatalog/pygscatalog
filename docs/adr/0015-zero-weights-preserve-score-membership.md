# Zero weights preserve score membership

## Status

Accepted.

## Decision

A matched scoring-file row belongs to its accession even when its effect weight
is zero. Its observed alleles contribute to `allele_count`, and its effect dosage
contributes to `dosage_sum`, exactly as for a nonzero weight. Its weighted score
contribution is zero. An all-zero score with called alleles therefore has sum zero
and average zero; an average with no observed alleles remains null.

This follows the biological definition of the denominator: it counts observed
alleles at matched scoring rows, not nonzero matrix coefficients. Adding a zero
coefficient does not make the genotype unobserved. A variant absent from an
accession is different from a variant explicitly present with a zero weight.

## Consequences

The dense weight matrix still uses zero for absent accession/variant pairs, which
is correct for multiplication. Statistics must derive membership from matched
scoring rows rather than testing the matrix coefficient for inequality to zero.
No new scoring option is needed.

The committed zero-weight v2.3.0 reference and an independent hard-call
calculation test this distinction, including an all-zero accession and a mixed
zero/nonzero accession. This corrects v3's previous omission of explicit zero
weights from its statistics; it does not copy v2's unrelated multi-score
denominator aggregation behaviour.
