# Extended scoring compatibility tests

These cases test selected alleles and calculated scores, not byte-for-byte output
or preservation of v2 defects. All normal test runs use committed inputs and
reference outputs without launching PLINK, Nextflow or Docker.

| Case | Agreement checked | Explicit difference or boundary |
| --- | --- | --- |
| Arithmetic | Multiple accessions, overlapping/disjoint variants, opposite effects, decimal signed weights, cancellation, additive/dominant/recessive contributions | v2 batches accessions into files with shared variant-union statistics; its final aggregate also includes denominators from unrelated effect-type files. v3 counts each accession's own observed alleles. Averages are independently checked under each denominator. |
| Partitioning | Reordered sample columns, reversed chromosome/score input order and complete sample batches reproduce the arithmetic reference | Fully called data only; changing an imputation population is not an invariance claim. MATCH receives all sample batches of a chromosome together. |
| Missingness | Called and missing GTs yield the same weighted sums under default v2 imputation; integer imputations 0, 1 and 2 depend on the called donors | v2 DENOM includes imputed calls while ALLELE_CT and NAMED_ALLELE_DOSAGE_SUM remain observed statistics. v3 averages use observed alleles and its dosage sum includes imputed values. Explicit v2 no-mean-imputation outputs document the policy difference. |
| Missingness rejection | An all-missing sample has a null v3 average; a wholly uncalled target is filtered then fails matching; a cohort below the imputation floor is rejected | These are specific v3 behaviours, not requirements to reproduce v2 failure wording. |
| Matching | Direct, swapped, complemented and no-other-allele matches; ambiguity and strand flags; exclusions, candidate selection and overlap rejection | Distinct ALT effects at one site and multi-base strand handling follow ADRs 0010/0011 rather than v2's ID collisions and plain complements. Selected target row and effect index are asserted explicitly. |
| Zero weights | An explicit zero-weight row retains its observed-allele and dosage statistics; an all-zero accession has a zero sum and zero average when calls exist | Corrects v3's previous use of nonzero coefficients as a membership mask; see ADR 0015. |

Counts and selected allele identities are exact. Decimal scores have an explicit
bound derived from PLINK's original six-significant-figure text precision, summed
over individual shards before cancellation. Independent GT/Decimal calculations
check v3 with a tighter floating-point bound; final averages allow its public
six-decimal rounding. Tolerances are not fitted to observed discrepancies.

Each reference contains original `.sscore` files, the final aggregate, matching
reports, software versions, run log, configuration and provenance. Tests check
the expected shard filenames, recorded input hashes and equivalent v2/v3 scoring
rows, including effect types. Generators reproduce inputs without writing
reference outputs. PAR, heteroplasmy, ancestry normalisation and fractional
dosages remain outside this suite; the original standalone PLINK case covers MT.

One further boundary emerged in the final source review: the current weight
matrix primary key contains target row and effect-allele index but not genetic
model. Different accessions using the same target/effect allele under different
models can therefore conflict. The mixed-model case here uses distinct scoring
rows; it does not establish support for that shared-allele combination.

## Deliberate regeneration

Use the unchanged upstream v2.3.0 commit and the recorded Docker image digests.
The script checks both before running. The work directory must be new and shared
with Docker; it must not be inside the pipeline checkout. The command writes the
selected case's reference directory, so review and commit its outputs together.

From the repository root:

```bash
git clone --branch v2.3.0 --depth 1 https://github.com/PGScatalog/pgsc_calc.git "$HOME/pgsc-calc-reference"
pgscatalog.utils/.venv/bin/python \
  pgscatalog.utils/packages/pgscatalog.calc/tests/data/v2_parity/generate_references.py \
  arithmetic --pipeline "$HOME/pgsc-calc-reference" \
  --work-dir "$HOME/pgsc-reference-arithmetic"
```

Run `arithmetic`, `missingness`, `zero_weight` with default mode. Also run
`missingness --mode no_mean_imputation`. Matching requires `matching` with each
of `--mode default`, `--mode keep_ambiguous` and `--mode ignore_strand_flips`.
Choose a distinct work directory for every invocation. `NXF_OFFLINE=true` is
optional when the required Nextflow components are cached.

The strand-exclusion mode uses the release's documented `task.ext.args`
configuration on MATCH_VARIANTS and MATCH_COMBINE. Merely passing an
`--ignore_strand_flips` Nextflow parameter does not wire the flag through this
release; the saved configuration and selected-match assertions verify the
effective behaviour. No workflow source is modified.

Run the suite with `nox -s v2_parity` from `pgscatalog.utils`, or use that
workspace's Python to run `pytest -c packages/pgscatalog.calc/pyproject.toml
packages/pgscatalog.calc/tests/integration/test_v2_parity*.py` from the same
directory.
