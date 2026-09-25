# Matching policy fixture

One indexed chromosome 1 VCF holds three fully called diploid samples (`S0`,
`S1`, `S2`) and 17 target records. Six additive scores isolate decisions that
would otherwise change the allele dosages before scoring. The extra score
position 999000 is absent from the VCF. All weights are exact small integers;
this fixture tests matching rather than floating-point precision.

| Accession | Question | v3 expectation |
| --- | --- | --- |
| PGS990020 | Direct REF, swapped ALT, complementary REF/ALT, missing other allele, incompatible alleles, absent position | Five matches; three when strand flips are disabled |
| PGS990021 | Two rows match the same REF allele through different other alleles | Exclude both collisions; retain the separate anchor |
| PGS990022 | Two rows match different ALT alleles of one multiallelic record | Retain both and the anchor; v2's ID-based collision rule differs |
| PGS990023 | Two target records at one position, with different genotypes | REFALT beats the earlier ALTREF candidate; retain the anchor |
| PGS990024 | Reverse complement, plain complement, inferred indel strand, forward indel | Match the reverse complement and forward indel; decline plain-complement and indel inference, per ADR 0011 |
| PGS990025 | Palindromic SNP and no-other-allele C at an A/C,G site | Exclude ambiguous pairings unless explicitly retained; the multiallelic ambiguity policy differs from v2 |

The anchor is an independent A/C SNP. It keeps accessions with deliberately
excluded rows scoreable at a 0.01 overlap threshold. A separate integration run
uses 0.8: only PGS990022 and PGS990023 survive, so exclusions cannot silently
become partial scores that bypass the acceptance gate.

Three v2 reference configurations use the same inputs: `default`,
`keep_ambiguous` and `ignore_strand_flips`, all with `--keep_multiallelic`.
The integration test compares shared-policy accessions with v2 and counts GT
alleles independently for intentional differences. It checks the actual target
record index and allele index used by scoring, not only log annotations. The
matching logs must account for all 23 scoring rows and name each exclusion.

The v2 workflow puts all six score columns in one PLINK scoring file. Its
`ALLELE_CT`, `DENOM` and `NAMED_ALLELE_DOSAGE_SUM` therefore describe the union
of selected alleles, including zero-filled weights for other accessions. v3
counts the retained rows of each score separately. Tests assert that specific
difference: v2 sums and its aggregate averages are checked against its union
denominator; v3 sums, counts, dosage sums and averages are checked against each
score's own selected GT alleles. Equal sums do not imply equal averages here.

`PGS*.txt` are v2 custom scoring files; `PGS*.txt.gz` carry identical rows in
v3's processed format. `generate.py` writes only these inputs and the indexed
VCF. Run it with the repository's Python environment to regenerate deterministic
inputs; it never writes reference outputs. The reproduction test compares all
14 generated files with the committed files byte for byte.
