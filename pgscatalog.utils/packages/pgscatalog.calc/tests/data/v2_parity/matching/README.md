# Matching policy fixture

One indexed chromosome 1 VCF contains three diploid samples: `S0`, `S1` and `S2`.
All samples have complete calls across 17 target records.
Six additive scores test matching decisions that can change allele dosages before scoring.
Scoring position 999000 is absent from the VCF.
All weights are small integers with exact values.
The fixture tests matching, not floating-point precision.

| Accession | Case | Expected v3 result |
| --- | --- | --- |
| PGS990020 | Direct REF, swapped ALT, complemented REF/ALT, missing other allele, incompatible alleles and absent position | Five matches. Three matches when strand flips are disabled. |
| PGS990021 | Two rows match the same REF allele through different other alleles | Exclude both colliding rows. Retain the separate anchor variant. |
| PGS990022 | Two rows match different ALT alleles of one multiallelic record | Retain both rows and the anchor variant. v2's ID-based collision rule differs. |
| PGS990023 | Two target records have the same position but different genotypes | Select REFALT instead of the earlier ALTREF candidate. Retain the anchor variant. |
| PGS990024 | Reverse complement, plain complement, inferred indel strand and forward indel | Match the reverse complement and forward indel. Reject the plain complement and indel strand inference, as specified in ADR 0011. |
| PGS990025 | A palindromic SNP and C without an other allele at an A/C,G site | Exclude ambiguous pairings unless the configuration retains them. The multiallelic ambiguity policy differs from v2. |

The anchor variant is a separate A/C SNP.
It lets accessions with deliberately excluded rows pass the 0.01 overlap threshold.
A separate integration run uses a threshold of 0.8.
Only PGS990022 and PGS990023 pass that threshold.
This check prevents partial scores from bypassing the acceptance gate after exclusions.

Three v2 reference configurations use the same inputs: `default`, `keep_ambiguous` and `ignore_strand_flips`.
All three use `--keep_multiallelic`.
The test compares accessions with shared policies against v2.
For intentional differences, it counts GT alleles independently.
It checks the target record index and allele index that scoring uses, in addition to the log annotations.
The matching logs must describe all 23 scoring rows and identify each exclusion.

## Statistics and averages

The v2 workflow places all six score columns in one PLINK scoring file.
Its `ALLELE_CT`, `DENOM` and `NAMED_ALLELE_DOSAGE_SUM` describe the union of selected alleles.
This union includes variants with zero-filled weights for other accessions.
v3 counts the retained rows of each score separately.

The tests check this difference explicitly:

- Check v2 sums and aggregate averages against its union denominator.
- Check v3 sums, counts, dosage sums and averages against each score's selected GT alleles.

Equal sums do not imply equal averages in this fixture.

## Regenerate inputs

`PGS*.txt` contains the v2 custom scoring files.
`PGS*.txt.gz` contains identical scoring rows in v3's processed format.
`generate.py` writes only these inputs and the indexed VCF.
Run it with the repository's Python environment to reproduce the inputs.
The generator never writes reference outputs.
The reproduction test compares all 14 generated files with the committed files byte for byte.
