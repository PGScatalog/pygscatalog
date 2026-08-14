# 0039 — `--sex_map`: the declaration, and when it is required

**Refs:** ADR 0014, ADR 0007, ADR 0003.

## Why

A **Declared sex** (`CONTEXT.md`) is the input the **Sex check** compares against, and it
is the half of ADR 0014 that needs no genotypes to test. Splitting it out keeps the
mandate, the parsing and the coverage rules provable on their own, and lands the fixture
churn — every chrX integration test gains a sex map — before the check exists to fail
them for a second reason.

Nothing here reads a genotype or decides anything about a sample. ADR 0007 still owns
ploidy; this ticket only carries a declaration to the door.

## Scope

**`--sex_map`, underscored**, in `add_load_input_args` (`cli/load.py:36`) alongside
`--bgen_sample_file` and `--absence_semantics`, so `pgsc_calc run` inherits it.

**Two columns, no header assumptions beyond sample ID and code.** Codes are `0` unknown,
`1` male, `2` female — plink and UK Biobank's convention. Any other value is refused,
naming the file, the line and the value. `0` is also how a user exempts a sample from the
check in 0040, so it is a legal value everywhere, never a parse warning.

**Mandatory when a scoring file carries sex chromosome variants.** Keyed on the scoring
files, not the target genomes: the Cache holds only positions the scoring files asked
for, so they are what decides whether there will be anything to check.
`Scorefiles.get_unique_positions(chrom=...)` (`lib/scorefile.py:70`) answers this, and
`Scorefiles` is already constructed at `cli/load.py:225` — move it above the genome
parsing so the mandate is enforced before any file is read. X and Y both trigger it.
**MT does not**: MT is hemizygous in everyone and says nothing about sex.

**Coverage is asymmetric, and the asymmetry is forced by sample batching (ADR 0003).**

- a genotyped sample **missing from the map** is a hard error — a mandatory declaration
  that is silently incomplete is not a declaration, and those samples would otherwise
  pass 0040 unchecked with nothing recorded
- a sample **in the map that this shard does not hold** is fine and expected. Batches see
  a subset of the cohort, so one cohort-wide map is the normal way to run this; rejecting
  extra rows would make the map a function of how the cohort was chopped up

Sample IDs come from the same place every other stage gets them, so a map keyed on the
wrong ID space fails as a coverage error rather than silently matching nothing.

**`pgsc_calc run` forwards it.** `cli/run.py:96` forwards `--bgen_sample_file`
conditionally and passes `--score_paths` into each LOAD dispatch, so without the same
three lines `run` cannot score any X-carrying PGS once the mandate lands.

**Existing fixtures gain a map**: `tests/integration/test_haploid_chrx.py`,
`tests/integration/test_pooled_imputation.py`, the chrX cases in
`tests/integration/test_bgen.py`, and every test in
`tests/integration/test_male_only_chry.py` — its `autosome_and_chry_scorefile` writes
chrY rows, and Y triggers the mandate as surely as X does, so those tests supply a map
and then take 0040's skip-and-warn path. This is the accepted cost of the mandate
(ADR 0014).

## Out of scope

The check itself — 0040.

**Reading sex from anywhere else.** No VCF header field records it. A bgen `.sample`
file's `sex` column is deliberately ignored even where one exists
(`tests/data/bgen/phased/tiny1000G.sample` has one): matching a column by name and
reading it under an assumed coding is two stacked inferences, and ADR 0014 keeps
inference out from under an exception that stops a run.

Deriving anything from the declaration. It never reaches a ploidy, a dosage or a score.

## Acceptance

- A scoring file carrying X or Y variants without `--sex_map` fails during argument
  validation, before any target genome is opened.
- A scoring file carrying MT variants and no X or Y does not require a map.
- Autosome-only scoring files are unaffected — the existing test suite passes untouched
  except for the chrX fixtures above.
- A map value outside `{0, 1, 2}` is refused naming the file, line and value.
- A genotyped sample absent from the map fails, naming the sample.
- A map covering the whole cohort loads cleanly against a shard holding one batch of it.
- `pgsc_calc run` scores an X-carrying PGS end to end, which fails without the forwarding.
