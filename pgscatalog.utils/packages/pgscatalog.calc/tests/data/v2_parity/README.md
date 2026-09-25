# PLINK 2 and pgsc_calc v2.3.0 parity fixture

The original all-chromosome case below remains unchanged. The additional
[`arithmetic`](arithmetic/README.md), [`missingness`](missingness/README.md),
[`matching`](matching/README.md) and `zero_weight` cases extend semantic coverage.
They use committed hard-call VCFs and original v2.3.0 outputs; fractional dosage
scenarios remain excluded. See [extended coverage](EXTENDED.md) for agreements,
intentional differences and reference regeneration.

This directory holds one self-contained fixture for comparing the public score
from the v3 scoring stages with two independent references: standalone PLINK 2
across all 25 chromosomes, and the actual pgsc_calc v2.3.0 workflow across its
24 supported chromosomes (1–22, X and Y). Both are configured to match v3's
scoring semantics. See [the v2.3.0 reference](pgsc_calc_v2.3.0/README.md) for its
coverage, configuration and regeneration instructions.

Add these files:

- `target.chr{chrom}.vcf.gz` and accompanying `.tbi` indexes: deterministic
  synthetic input split into 25 VCFs, one per chromosome: `1`–`22`, `X`, `Y`,
  and `MT`. Use canonical chromosome names inside the VCFs.
- `score.txt.gz`: the harmonized scoring file used by v3.
- `plink.weights.tsv`: equivalent variant IDs, effect alleles, and weights
  for PLINK 2.
- `sex.tsv`: v3 sample-sex map.
- `samples.psam` and `samples.Y.psam`: PLINK sample metadata for the full
  cohort and the XY-only Y file, respectively.
- `plink_reference/chr{chrom}.sscore`: the original standalone PLINK 2 results.
- `pgsc_calc_v2.3.0/`: raw chromosome scores, final aggregate and provenance
  from the actual v2.3.0 workflow.

The generator produces the VCFs, indexes, scoring inputs, and sample-sex
files from one fixed variant and sample definition. Running PLINK 2 and
refreshing the saved reference are separate, explicit steps.
The generator must not overwrite either reference directory. The integration
test sums committed chromosome score sums by sample ID using Decimal and rounds
only the final totals to six decimal places. It also checks v2.3.0's original
aggregate against that workflow's chromosome outputs.

The fixture contains four samples in fixed order: `XX_1`, `XX_2`, `XY_1`,
`XY_2`. Use explicit, fixed genotypes that collectively exercise effect-allele
dosages 0, 1, and 2. Include an ALT effect allele, a REF effect
allele, a positive non-integer weight, and a negative non-integer weight.
Include biallelic and multiallelic variants across every human autosome and
both allosomes (X and Y), plus MT. Each chromosome contains a biallelic SNP
and a multiallelic SNP with two ALT alleles. Chromosome 1 additionally contains
one short insertion and one short deletion (52 records in total). Use identical
indel positions and allele strings in the VCF and scoring inputs; testing
normalization of equivalent representations is outside this fixture's scope.
Allele sequences are synthetic and contain only uppercase `A`, `C`, `T`, and
`G`. Indels use explicit nucleotide sequences with a shared anchor base; do
not use `N`, `*`, symbolic alleles, or I/D codes. Positions and alleles are
fixed, with no reference FASTA dependency or claim that REF matches a real
reference sequence. Reference-sequence validation and liftover are outside
this fixture's scope.
Use one additive score spanning all 25 chromosomes, with one scored allele
per variant. Alternate REF and ALT effect alleles across records, selecting
the second ALT at some multiallelic records. Use fixed positive and negative
fractional weights. The expected output has four rows, one per sample.
Fixed calls include `1/2` at diploid
multiallelic variants and individual ALT alleles at haploid variants.
Every genotype is fully called; missing genotypes and imputation
are outside this fixture's scope. MT genotypes are haploid, with no
heteroplasmy.

Exclude pseudoautosomal regions (PAR) from this fixture. All X variants are
non-pseudoautosomal, with diploid XX calls and haploid XY calls. Testing
ploidy changes across PAR boundaries is deferred. Use GRCh38 as the coordinate
convention, without claiming reference-sequence agreement; choose X and Y
positions outside its PARs.

All autosomal, X, and MT files share the same samples. The Y file contains only
XY samples, with fully called haploid genotypes: the current VCF reader
rejects missing Y calls because it cannot distinguish an absent chromosome
from a missing genotype. Its sample columns are `XY_1`, `XY_2`, in that order.
Supply explicit sex metadata for both tools: `XX_1` and `XX_2` have code 2;
`XY_1` and `XY_2` have code 1. Chromosome files belong to one sampleset, and their score
contributions are combined across chromosomes.

Configure PLINK 2 explicitly with `--xchr-model 1` so a male non-pseudoautosomal
X effect allele contributes dosage 1, matching v3. Record any additional
import and scoring options required to align the two runs.

## Generate and check

Run from `pgscatalog.utils` with its development environment installed:

```bash
fixture=packages/pgscatalog.calc/tests/data/v2_parity
uv run python "$fixture/generate.py" "$fixture"
uv run pytest -q -c packages/pgscatalog.calc/pyproject.toml "$fixture/test_generate.py"
```

`generate.py` uses the standard library and the existing `pysam` dependency.
It overwrites only its 55 generated inputs, never `plink_reference/`. There is no
random seed: every allele, genotype, position, ID and weight is fixed. Gzip
headers have no filename or timestamp. Regeneration is byte-identical within
the same Python/pysam/compression-library environment. The check regenerates
in two directories, compares bytes, fetches through every index, and checks
sample membership, allele counts, ploidy and scoring-input consistency.

Weights are multiples of 1/16, exactly representable in binary floating point.
Their magnitude varies by chromosome so missing chromosome contributions are
easier to detect. The accession is the synthetic `PGS999997`, not a real score.

## Produce the PLINK reference

Run this Bash block from the same directory, with `plink2` on PATH. Use a fresh
output directory and retain its logs, version and unmodified `.sscore` files.

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

`--split-par b38` satisfies PLINK's explicit X-import handling; it changes no
records because the fixture has no PAR variants. Its no-effect warning is
expected. Each chromosome scores only its own IDs from the shared weights
file, so warnings about other chromosomes' missing IDs are also expected.
There must be four scored variants on chromosome 1 and two on every other
chromosome. Do not use `WEIGHT_AVG` as the score: combine `WEIGHT_SUM`.

After checking the PLINK logs, copy the 25 original `.sscore` files into
`plink_reference/` and refresh the provenance and saved logs. Confirm all 52
variant IDs were scored using the `.sscore.vars` files before recording them
in `provenance.json`. The integration test performs the aggregation itself;
there is no separate aggregation command or generated expected-score CSV.

## Run the current v3 stages

`pgsc_calc run` does not currently expose `--keep_multiallelic`; use the stage
commands below. The current match-rate gate uses strict `>` rather than `>=`,
so `min_overlap=1` rejects even complete matching. Until that boundary is
fixed, use `0.99` for this 52-variant fixture: 51/52 is below 0.99, so only all
52 matching variants can pass. The intended integration-test setting remains
`1` once the boundary is fixed.

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

The result is `$v3_out/combined/scores.txt.gz`. The two X variants exercise
scoring ploidy, but are too few for a conclusive sex check; that warning is
expected. PAR handling, missingness, imputation, strand-flip cases, ambiguous
matches, dominant/recessive scores and reference normalisation are deferred.

## Saved reference and integration test

The committed `.sscore` files were generated with **PLINK v2.0.0-a.6.23 M1 (24 Aug 2025)**.
The inputs were regenerated with Python 3.12.8 and pysam 0.23.3 using
`python generate.py .` in this directory. The reference uses the commands
above, with temporary output prefixes. The exact argument lists, versions,
input SHA-256 hashes and 52 scored variant IDs are saved in
[`plink_reference/provenance.json`](plink_reference/provenance.json).

[`plink_reference/`](plink_reference/) retains the 25 original chromosome
`.sscore` outputs and the concatenated original PLINK logs in `plink.log`.
The temporary output paths in those logs are historical; the corresponding
score files are retained here under the same basenames. The integration test sums
`WEIGHT_SUM` by `#IID` with Decimal and rounds only the final sum to six decimal
places. It validates chromosome file coverage and sample membership before
aggregation, rejects duplicate sample rows and non-finite contributions. XX samples have no Y
contribution. Neither v3 nor the generator produced these PLINK results.

The test also sums PLINK's `ALLELE_CT` and `NAMED_ALLELE_DOSAGE_SUM` to check
v3's `allele_count` and `dosage_sum`. It derives `score_avg` by dividing the
combined `WEIGHT_SUM` by the combined `ALLELE_CT`, then rounds once to six
decimal places; chromosome averages are never averaged. Each chromosome's
partial allele counts are checked separately against PLINK because the X/Y
contributions happen to give all four samples the same total of 98 alleles,
which could otherwise hide errors that cancel in the combined denominator.

The integration test loads all 25 checked-in VCFs, supplies the sex map, keeps
multiallelic matches, and runs the public v3 stages through COMBINE. It requires
all 52 variants to match, all four sample/accession keys to be present without
duplicates, and exact equality of the public six-decimal score sums against the totals it
calculates directly from the committed PLINK outputs. It uses
the documented `min_overlap=0.99` workaround. It never invokes PLINK or Nextflow,
generates fixtures, or rewrites the reference.

Run only the integration test from `pgscatalog.utils`:

```bash
nox -s v2_parity
```

Or run both the integration test and generator check in the existing environment:

```bash
uv run pytest -q -c packages/pgscatalog.calc/pyproject.toml \
    packages/pgscatalog.calc/tests/integration/test_v2_parity.py \
    packages/pgscatalog.calc/tests/data/v2_parity/test_generate.py
```

The integration test has two cases: standalone PLINK 2 (52 variants) and
pgsc_calc v2.3.0 (50 variants; MT remains covered by the first case).
Only refresh the raw `.sscore` files and their provenance after inspecting and understanding
any change to the results.
