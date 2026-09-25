# PLINK 2 and pgsc_calc v2.3.0 parity fixture

This fixture compares scores from the public v3 stages with two external references.
Both references use PLINK. The standalone PLINK 2 run covers 25 chromosomes.
The pgsc_calc v2.3.0 workflow covers its 24 supported chromosomes: 1–22, X and Y.
Both runs use settings that match v3's scoring policy for this fixture.
See the [v2.3.0 reference](pgsc_calc_v2.3.0/README.md) for its settings and regeneration procedure.

The original fixture remains unchanged. Additional cases cover
[arithmetic](arithmetic/README.md), [missingness](missingness/README.md),
[matching](matching/README.md) and [zero weights](zero_weight/README.md).
These cases use committed VCF hard calls and original v2.3.0 outputs.
They exclude fractional dosage scenarios.
See [extended coverage](EXTENDED.md) for agreements, intentional differences and regeneration instructions.

## Files

- `target.chr{chrom}.vcf.gz`: 25 synthetic VCFs, one for each chromosome from `1`–`22`, `X`, `Y` and `MT`.
  Each VCF uses canonical chromosome names and has a `.tbi` index.
- `score.txt.gz`: the harmonised scoring file for v3.
- `plink.weights.tsv`: equivalent variant IDs, effect alleles and weights for PLINK 2.
- `sex.tsv`: the sample-sex map for v3.
- `samples.psam`: PLINK metadata for all four samples.
- `samples.Y.psam`: PLINK metadata for the two XY samples in the Y file.
- `plink_reference/chr{chrom}.sscore`: the original standalone PLINK 2 results.
- `pgsc_calc_v2.3.0/`: the original chromosome scores, final aggregate and provenance from the v2.3.0 workflow.

The generator uses one fixed definition of the variants and samples.
It creates the VCFs, indexes, scoring inputs and files that declare sample sex.
It must not overwrite either reference directory.
Reference generation requires a separate, deliberate run.

The integration test sums the committed chromosome scores by sample ID with `Decimal`.
It rounds the final totals to six decimal places.
It also checks the original v2.3.0 aggregate against that workflow's chromosome outputs.

## Variants and samples

The sample order is `XX_1`, `XX_2`, `XY_1`, `XY_2`.
Fixed genotypes provide effect-allele dosages of 0, 1 and 2.
Each chromosome has a biallelic SNP and a multiallelic SNP with two ALT alleles.
Chromosome 1 also has one short insertion and one short deletion.
The fixture contains 52 records across all autosomes, both allosomes (X and Y), and MT.

Allele sequences contain only uppercase `A`, `C`, `T` and `G`.
Each indel uses explicit nucleotide sequences with a shared anchor base.
The fixture excludes `N`, `*`, symbolic alleles and I/D codes.
The VCFs and scoring inputs use identical indel positions and allele strings.
The fixture does not test normalisation of equivalent indel representations.

Positions and alleles are fixed. The fixture has no reference FASTA dependency.
It does not claim that REF matches a real reference sequence.
Reference-sequence checks and liftover remain outside its scope.

The fixture uses one additive score, `PGS999997`, with one scored allele per variant.
This accession identifies a synthetic score.
The effect allele alternates between REF and ALT across records.
Some multiallelic records use the second ALT as the effect allele.
Weights include fixed positive and negative non-integer values.
The expected output has one row per sample, for a total of four rows.

Fixed diploid calls include `1/2` at multiallelic variants.
Haploid calls include individual ALT alleles.
Every genotype has a complete call.
Missing genotypes and imputation remain outside this original fixture.
MT genotypes are haploid and have no heteroplasmy.

## Ploidy and sex metadata

The fixture excludes pseudoautosomal regions (PAR).
It uses GRCh38 coordinates, with all X and Y positions outside the PARs.
X calls are diploid in XX samples and haploid in XY samples.
The fixture does not test ploidy changes across PAR boundaries.

Autosomal, X and MT files contain all four samples.
The Y file contains only `XY_1` and `XY_2`, in that order, with complete haploid calls.
The current VCF reader rejects missing Y calls.
It cannot distinguish an absent chromosome from a missing genotype.

Both tools receive explicit sex metadata.
`XX_1` and `XX_2` have code 2. `XY_1` and `XY_2` have code 1.
All chromosome files belong to one sampleset.
The test combines their score contributions across chromosomes.

Configure PLINK 2 with `--xchr-model 1`.
This setting gives a male non-PAR X effect allele dosage 1, which matches v3.
Record any other import or scoring options required to align the runs.

## Generate and check inputs

Install the development environment before you run these commands.
Run them from `pgscatalog.utils`:

```bash
fixture=packages/pgscatalog.calc/tests/data/v2_parity
uv run python "$fixture/generate.py" "$fixture"
uv run pytest -q -c packages/pgscatalog.calc/pyproject.toml "$fixture/test_generate.py"
```

`generate.py` uses the standard library and the existing `pysam` dependency.
It overwrites only its 55 generated inputs.
It never writes to a reference directory.
Every allele, genotype, position, ID and weight is fixed, so the generator needs no random seed.
Gzip headers contain no filename or timestamp.
Regeneration produces identical bytes within the same Python, pysam and compression-library environment.

The generator test performs these checks:

1. Regenerate inputs in two directories.
2. Compare the generated bytes.
3. Read records through every index.
4. Check sample membership, allele counts and ploidy.
5. Check agreement between the scoring inputs.

Weights are multiples of 1/16, which binary floating-point values can represent exactly.
Their magnitudes vary by chromosome. This variation helps the test detect a missing chromosome contribution.

## Produce the PLINK reference

Run this Bash block from `pgscatalog.utils` with `plink2` on PATH.
Use a new output directory.
Retain the logs, version and original `.sscore` files.

```bash
set -euo pipefail
fixture=packages/pgscatalog.calc/tests/data/v2_parity
plink_out=$(mktemp -d)
plink2 --version > "$plink_out/version.txt"
for chrom in {1..22} X Y MT; do
    metadata="$fixture/samples.psam"
    extra=()
    if [[ "$chrom" == Y ]]; then metadata="$fixture/samples.Y.psam"; fi
    if [[ "$chrom" == X ]]; then extra=(--split-par b38); fi
    plink2 --vcf "$fixture/target.chr$chrom.vcf.gz" \
        --psam "$metadata" "${extra[@]}" --xchr-model 1 \
        --score "$fixture/plink.weights.tsv" 1 2 3 header-read \
        no-mean-imputation list-variants cols=+scoresums \
        --threads 1 --out "$plink_out/chr$chrom"
done
```

`--split-par b38` satisfies PLINK's requirement for explicit X import handling.
It changes no records because the fixture has no PAR variants.
The warning that it has no effect is expected.
Each chromosome uses only its own IDs from the shared weights file.
Warnings about missing IDs from other chromosomes are also expected.

Chromosome 1 must have four scored variants. Every other chromosome must have two.
Combine `WEIGHT_SUM` values to calculate the total score.
Do not use `WEIGHT_AVG` as the score sum.

1. Check the PLINK logs.
2. Check the `.sscore.vars` files for all 52 scored variant IDs.
3. Copy the 25 original `.sscore` files into `plink_reference/`.
4. Update `provenance.json` with the scored IDs and details of the run.
5. Update the saved logs.

The integration test calculates the aggregate itself.
It needs no separate aggregation command or generated CSV of expected scores.

## Run the current v3 stages

`pgsc_calc run` does not currently expose `--keep_multiallelic`.
Use the separate stage commands below.

The current match-rate gate uses `>` instead of `>=`.
Thus, `min_overlap=1` rejects even a complete match.
Use `0.99` for this 52-variant fixture until the gate accepts equality.
A result of 51/52 is below 0.99, so only 52 matches pass.
The intended test setting is `1` after the boundary correction.

```bash
set -euo pipefail
fixture=packages/pgscatalog.calc/tests/data/v2_parity
v3_out=$(mktemp -d)
for chrom in {1..22} X Y MT; do
    uv run pgsc_calc load \
        --target_genomes "$fixture/target.chr$chrom.vcf.gz:$chrom" \
        --sampleset synthetic --format vcf --score_paths "$fixture/score.txt.gz" \
        --sex_map "$fixture/sex.tsv" --cache_dir "$v3_out/load/$chrom"
    uv run pgsc_calc match \
        --zarr_zip_file "$v3_out/load/$chrom/genotypes.zarr.zip" \
        --score_paths "$fixture/score.txt.gz" --keep_multiallelic \
        --out_dir "$v3_out/match/$chrom"
done
uv run pgsc_calc match-combine \
    --matches "$v3_out"/match/*/matches.parquet \
    --samples "$v3_out"/match/*/samples.parquet \
    --score_paths "$fixture/score.txt.gz" --sampleset synthetic \
    --min_overlap 0.99 --out_dir "$v3_out/matched"
for chrom in {1..22} X Y MT; do
    uv run pgsc_calc score \
        --zarr_zip_file "$v3_out/load/$chrom/genotypes.zarr.zip" \
        --score_paths "$fixture/score.txt.gz" \
        --match_table "$v3_out/matched/match_table.parquet" \
        --out_dir "$v3_out/score/$chrom"
done
uv run pgsc_calc combine --partial_scores "$v3_out"/score/*/partial_scores.parquet \
    --out_dir "$v3_out/combined"
```

The result is `$v3_out/combined/scores.txt.gz`.
The two X variants test scoring ploidy, but they are insufficient for a conclusive sex check.
The warning about this limit is expected.

The original fixture excludes PAR handling, missingness, imputation, strand flips, ambiguous matches, dominant/recessive scores and reference normalisation.
The [extended cases](EXTENDED.md) test some of these behaviours separately.

## Saved reference and assertions

The standalone reference uses **PLINK v2.0.0-a.6.23 M1 (24 Aug 2025)**.
The inputs were regenerated with Python 3.12.8 and pysam 0.23.3, using `python generate.py .` in this directory.
The reference uses the commands above with temporary output prefixes.
[`plink_reference/provenance.json`](plink_reference/provenance.json) records the arguments, versions, input SHA-256 hashes and 52 scored IDs.

[`plink_reference/`](plink_reference/) contains the 25 original `.sscore` files.
`plink.log` contains the original PLINK logs in sequence.
The temporary paths in those logs identify the original run.
The saved score files retain their original basenames.
Neither v3 nor the generator produced these reference scores.

The test checks the expected chromosome files and samples before it combines scores.
It rejects duplicate sample rows and non-finite contributions.
XX samples have no Y contribution.

The test uses these calculations:

- Sum `WEIGHT_SUM` by `#IID` with `Decimal`.
- Sum `ALLELE_CT` to check v3's `allele_count`.
- Sum `NAMED_ALLELE_DOSAGE_SUM` to check v3's `dosage_sum`.
- Divide the combined `WEIGHT_SUM` by the combined `ALLELE_CT` to calculate `score_avg`.
- Round each final score and average once to six decimal places.

The test does not average chromosome averages.
It checks each chromosome's score, dosage sum and allele count before it checks the totals.
All four samples have 98 alleles in the full fixture because the X and Y counts balance.
Separate chromosome checks detect errors that could cancel in a total.

The standalone case loads all 25 committed VCFs and supplies the sex map.
It retains multiallelic matches and runs the public v3 stages through COMBINE.
It requires all 52 variants to match.
It checks four unique sample/accession keys and exact equality of the final six-decimal score sums.
It uses the `min_overlap=0.99` workaround described above.
The test does not launch PLINK or Nextflow, generate inputs, or rewrite references.

Run the full parity suite from `pgscatalog.utils`:

```bash
nox -s v2_parity
```

Alternatively, run the original integration test and generator check in the existing environment:

```bash
uv run pytest -q -c packages/pgscatalog.calc/pyproject.toml \
    packages/pgscatalog.calc/tests/integration/test_v2_parity.py \
    packages/pgscatalog.calc/tests/data/v2_parity/test_generate.py
```

The original integration test has two cases.
The standalone PLINK 2 case covers 52 variants.
The pgsc_calc v2.3.0 case covers 50 variants, with MT covered by the standalone case.
Inspect and understand changes to results before you replace raw `.sscore` files or their provenance.
