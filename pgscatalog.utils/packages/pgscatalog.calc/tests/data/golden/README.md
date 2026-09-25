# Golden scores

The original `scores.txt.gz` contains results from one process that scored the fixture data.
That run preceded the split into MATCH, MATCH_COMBINE, SCORE and SCORE_COMBINE described in ADR 0003.
SCORE now receives only one shard, so this file preserves the earlier reference.

Every target layout in `tests/integration/conftest.py` contains the same sample/variant pairs.
Each layout must therefore reproduce these numbers.

## Correction for zero weights

ADR 0015 corrected the statistics on 25 September 2026.
The fixture contains 123 matched diploid scoring rows.
Six rows have zero weights: row numbers 10, 30, 50, 70, 90 and 110.
The earlier calculation excluded these rows from allele counts and dosage statistics.

All 3,202 samples have complete calls.
Including the six rows changes each allele count from 234 to 246.
Dosage sums now include observed effect-allele copies at these six loci.
Averages now use 246 as the denominator.
Weighted sums and `n_matched` remain unchanged because a dosage multiplied by zero contributes nothing to the score.

The correction used `unsplit_target` VCF GT calls and `scorefile` fixture rows.
It did not read v3's calculated output.
Chromosome, position, REF and ALT together identified each variant.
Position alone is not unique in this fixture.

The independent calculation used these checks:

1. Count effect-allele copies directly from each called diploid GT for each sample.
2. Calculate each weighted sum with `Decimal`.
3. Check that the sum equals the original golden sum before writing an output.
4. Check that the nonzero-weight rows reproduce every original allele count and dosage sum.

The corrected file preserves the original score strings and match counts.
It rounds corrected averages to six decimal places with `ROUND_HALF_UP`.
The gzip header uses a zero timestamp and no filename to keep the output deterministic.

## Regeneration

Run `pgsc_calc run` on the `unsplit_target` fixture.
See `run_serial` in `tests/integration/conftest.py` for the invocation.
Inspect and understand any intended changes to the results.
Then copy the resulting `scores.txt.gz` into this directory.
