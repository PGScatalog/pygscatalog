# Zero weights preserve score membership

## Status

Accepted.

## Decision

A matched scoring-file row belongs to its accession even when its effect weight
is zero. Its observed alleles contribute to `allele_count`. Its effect dosage
contributes to `dosage_sum`, exactly as for a nonzero weight. Its weighted score
contribution is zero.

An all-zero score with called alleles therefore has sum zero and average zero.
An average with no observed alleles remains null.

The denominator counts observed alleles at matched scoring rows, not nonzero
matrix coefficients. A zero coefficient does not make the genotype unobserved.
A variant absent from an accession differs from a variant explicitly present
with a zero weight.

## Consequences

The dense weight matrix still uses zero for absent accession/variant pairs. This
is correct for multiplication. Statistics must derive membership from matched
scoring rows, not from whether the matrix coefficient differs from zero. No new
scoring option is needed.

The committed v2.3.0 zero-weight reference and an independent hard-call
calculation test this distinction. They include an all-zero accession and a mixed
zero/nonzero accession. This corrects v3's previous omission of explicit zero
weights from its statistics. It does not copy v2's unrelated behaviour when
aggregating denominators across multiple scores.
