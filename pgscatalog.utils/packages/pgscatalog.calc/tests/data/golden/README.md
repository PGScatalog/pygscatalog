# Golden scores

`scores.txt.gz` is the output of the fixture data scored by a single process, before
SCORE was split into MATCH / MATCH_COMBINE / SCORE / SCORE_COMBINE (docs/adr/0003).
Once `score` only ever sees one shard there is no single-process reference left, so
this file is it.

Every target genome shape in `tests/integration/conftest.py` holds the same
(sample, variant) pairs, so all of them must reproduce these numbers.

## Zero-weight membership correction

ADR 0015 corrected the statistics on 25 September 2026. The fixture has 123
matched diploid scoring rows, six with zero effect weights (row numbers 10, 30,
50, 70, 90 and 110). The earlier calculation omitted these six rows from allele
counts and dosage statistics. All 3,202 samples have complete calls, so including
them changes each allele count from 234 to 246. Dosage sums now include the
observed effect-allele copies at these six loci, and averages divide by 246.
Weighted sums and `n_matched` are unchanged: multiplying a dosage by zero cannot
change the score.

This migration was calculated independently from the `unsplit_target` VCF GT
calls and `scorefile` fixture rows, without reading v3's calculated output.
Variants were identified by chromosome, position, REF and ALT (position alone
is not unique in this fixture). For each sample, the effect allele's copies
were counted directly from each called diploid GT; their Decimal-weighted sum
was checked against the original golden sum before any output was written.
Restricting to nonzero weights also reproduced every original allele count and
dosage sum. The corrected file preserves the original score strings and match
counts, and rounds the corrected averages to six decimals using ROUND_HALF_UP.
The gzip header uses a zero timestamp and no filename for deterministic output.

To regenerate, run `pgsc_calc run` on the `unsplit_target` fixture (see
`run_serial` in `tests/integration/conftest.py`) and copy the resulting
`scores.txt.gz` here. Only do this if the change to the numbers is intended and
understood.
