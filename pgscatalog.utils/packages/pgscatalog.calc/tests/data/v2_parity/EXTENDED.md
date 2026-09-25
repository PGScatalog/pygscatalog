# Extended scoring compatibility tests

These cases test selected alleles and calculated scores.
They do not require identical output bytes or preservation of v2 defects.
Normal tests read committed inputs and reference outputs.
They do not launch PLINK, Nextflow or Docker.

## Coverage and differences

### Arithmetic

The tests cover multiple accessions with overlapping and disjoint variants, opposite effects, decimal weights and cancellation.
They also cover additive, dominant and recessive contributions.

v2 groups accessions into files with shared statistics for the union of their variants.
Its final aggregate also includes denominators from unrelated files for other effect types.
v3 counts the observed alleles for each accession separately.
The tests check each version's averages with that version's denominator.

### Partitioning

The tests change sample-column order, chromosome-input order and scoring-file order.
They also divide the samples into complete batches.
All results must agree with the arithmetic reference.
MATCH receives all sample batches for a chromosome together.

These tests use complete genotype calls.
They do not claim that results remain unchanged when a change to sample batches changes the population used for imputation.

### Missingness

Default v2 imputation and v3 produce the same weighted sums for the fixture's called and missing genotypes.
Called genotypes determine integer imputed dosages of 0, 1 and 2.

v2's DENOM includes imputed calls.
Its ALLELE_CT and NAMED_ALLELE_DOSAGE_SUM describe only observed calls.
v3 uses observed alleles as the average denominator, but its dosage sum includes imputed values.
A separate v2 reference disables mean imputation to show the difference between these policies.

The tests also check specific v3 outcomes:

- A sample with no called genotypes has a null average.
- LOAD removes a wholly uncalled target, which then fails matching.
- v3 rejects a cohort below the minimum sample count for imputation.

The tests do not require identical error wording in v2 and v3.

### Matching

The tests cover direct, swapped and complemented matches, including matches without an other allele.
They check ambiguity flags, strand flags, exclusions, candidate selection and overlap rejection.
They also check the selected target row and effect-allele index.

v3 follows ADRs 0010 and 0011 for distinct ALT effects at one site and strand handling for multi-base alleles.
It does not reproduce v2's ID collisions or plain-complement behaviour where those policies differ.

### Zero weights

An explicit zero-weight row retains its observed-allele and dosage statistics.
A score with only zero weights has a zero sum and zero average when calls exist.
ADR 0015 corrects v3's previous use of nonzero coefficients to identify score membership.

## Numerical checks

The tests require exact counts and selected allele identities.
PLINK writes the reference scores with six significant figures.
The permitted error follows that text precision.
The tests add the error bounds for individual shards before score contributions cancel.
A shard is one part of the data used for scoring.

Independent calculations from GT calls and `Decimal` weights check v3 with a smaller floating-point error bound.
Checks of final averages allow v3's rounding to six decimal places.
The tests do not choose tolerances from observed differences.

## Reference files and limits

Each reference contains these files:

- Original `.sscore` outputs and the final aggregate.
- Matching reports and software versions.
- The run log, configuration and provenance.

Tests check the expected shard filenames and recorded input hashes.
They check equivalent v2 and v3 scoring rows, including effect types.
Generators reproduce inputs without writing reference outputs.

The suite excludes PAR, heteroplasmy, ancestry normalisation and fractional dosages.
The original standalone PLINK case covers MT.

The current weight-matrix primary key contains the target row and effect-allele index, but not the genetic model.
Different accessions can therefore conflict if they use the same target/effect allele with different models.
The mixed-model case uses distinct scoring rows.
It does not establish support for that shared-allele combination.

## Regenerate references

Use the unchanged upstream v2.3.0 commit and recorded Docker image digests.
The script checks both before it runs the workflow.
Use a new work directory that Docker can read.
Place that directory outside the pipeline checkout.
The command writes to the selected case's reference directory.
Review the outputs together before you commit them.

Run from the repository root:

```bash
git clone --branch v2.3.0 --depth 1 https://github.com/PGScatalog/pgsc_calc.git "$HOME/pgsc-calc-reference"
pgscatalog.utils/.venv/bin/python \
  pgscatalog.utils/packages/pgscatalog.calc/tests/data/v2_parity/generate_references.py \
  arithmetic --pipeline "$HOME/pgsc-calc-reference" \
  --work-dir "$HOME/pgsc-reference-arithmetic"
```

Generate each reference with a separate work directory:

| Case | Mode |
| --- | --- |
| `arithmetic` | Default |
| `missingness` | Default |
| `missingness` | `--mode no_mean_imputation` |
| `zero_weight` | Default |
| `matching` | `--mode default` |
| `matching` | `--mode keep_ambiguous` |
| `matching` | `--mode ignore_strand_flips` |

`NXF_OFFLINE=true` is optional when the required Nextflow components are cached.

The strand-exclusion mode uses `task.ext.args` on MATCH_VARIANTS and MATCH_COMBINE.
This release does not pass a top-level `--ignore_strand_flips` Nextflow parameter to those commands.
The saved configuration and selected-match assertions check the effective behaviour.
The workflow source remains unchanged.

Run `nox -s v2_parity` from `pgscatalog.utils` to test the suite.
Alternatively, use the workspace Python from that directory with this command:

`pytest -c packages/pgscatalog.calc/pyproject.toml packages/pgscatalog.calc/tests/integration/test_v2_parity*.py`
