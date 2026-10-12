# pgsc_calc v2.3.0 reference

The upstream Nextflow workflow produced these outputs at tag `v2.3.0`, commit `72ee54f2119d2fe65ac671327342dfc412bd0a3b`.
The run used the release's pinned Docker images and unchanged workflow source.
`provenance.json` records image digests, the command, original input hashes and chromosome scoring scripts.
`software_versions.yml` contains the workflow's version report.

## Coverage and configuration

The pinned matcher, `pgscatalog.match` 0.4.0, accepts target chromosomes 1–22, X and Y.
This reference therefore covers **24 chromosomes and 50 variants**.
Its samplesheet excludes MT.
`PGS999997.txt` contains all 52 original weights, positions and allele pairs in v2's custom scoring format.
The matching summary and full log record two unmatched MT rows.
`--min_overlap 0.95` accepts 50/52 matches and rejects 49/52 matches.
The integration test also requires exactly 50 matched variants.

The test compares selected matches by accession and scoring-file row.
Each match must have the same chromosome, position and effect allele in v2 and v3.
The test also checks v3's effect-allele indices against the original VCF alleles.
It excludes unselected v2 candidates, duplicate flags and ALT annotations from this comparison.
v2 joins records by ID alone. This join can duplicate multiallelic candidates and attach an incorrect ALT annotation.
The effect allele sent to PLINK can remain correct despite this annotation error.

VCF import receives the committed sample-sex metadata.
Scoring uses `--xchr-model 1` to match v3 and the standalone PLINK reference.
`--split-par b38` changes no variants because all X positions are outside PAR.
`--keep_multiallelic` retains multiallelic variants.
The run requests uncompressed `.sscore` files to simplify review.
These results describe this configuration. They do not establish equivalence with v2 defaults.

`--only_score` excludes ancestry calculation and report generation.
The configuration requests reduced resources and one CPU for the small input.
The run used Nextflow 25.10.0 and the arm64 images recorded in `provenance.json`.

## Committed outputs

- `synthetic_{chrom}_additive_0.sscore`: 24 unchanged PLINK outputs from the v2 workflow.
  The test sums `PGS999997_SUM`, `ALLELE_CT` and `NAMED_ALLELE_DOSAGE_SUM` by `IID`.
  It then calculates the score average.
- `aggregated_scores.txt.gz`: the original SCORE_AGGREGATE output.
  The test checks its SUM, DENOM and AVG against the chromosome results.
- `synthetic_summary.csv` and `synthetic_log.csv.gz`: the original matching reports.
- `plink.log`: the original scoring logs in sequence.
- `nextflow.log`: the workflow console output, which records successful completion.
- `PGS999997.txt`, `samplesheet.csv` and `reference.config`: the original inputs.
  The samplesheet and configuration retain the original absolute paths.

The v3 comparison uses the same committed scoring data and 24 VCFs.
The separate standalone PLINK case covers all 25 chromosomes, including MT.
Neither test launches PLINK, Nextflow or Docker.

## Regeneration

Nextflow and Docker must be available.
Run these commands from `pgscatalog.utils`.
They create a new directory and update only the saved absolute fixture path.

```bash
fixture=$(pwd)/packages/pgscatalog.calc/tests/data/v2_parity
v2_run=$(mktemp -d "$HOME/pgsc-calc-v230.XXXXXX")
python - "$fixture" "$v2_run" <<'PY'
import csv
import shutil
import sys
from pathlib import Path

fixture, output = map(Path, sys.argv[1:])
source = fixture / "pgsc_calc_v2.3.0"
with (source / "samplesheet.csv").open() as handle:
    original = Path(next(csv.DictReader(handle))["path_prefix"]).parent
for name in ("samplesheet.csv", "reference.config"):
    text = (source / name).read_text().replace(str(original), str(fixture))
    (output / name).write_text(text)
shutil.copyfile(source / "PGS999997.txt", output / "PGS999997.txt")
PY
git clone --depth 1 --branch v2.3.0 https://github.com/PGScatalog/pgsc_calc.git "$v2_run/pipeline"
test "$(git -C "$v2_run/pipeline" rev-parse HEAD)" = 72ee54f2119d2fe65ac671327342dfc412bd0a3b
cd "$v2_run"
nextflow run pipeline -profile docker -c reference.config \
    --input samplesheet.csv --scorefile PGS999997.txt --target_build GRCh38 \
    --keep_multiallelic --min_overlap 0.95 --only_score --outdir results \
    -ansi-log false > run.log 2>&1
```

Use a directory that the Docker VM can read.
On this host, Docker could not read `/private/tmp`.
The recorded run used `NXF_OFFLINE=true` because the Nextflow plugin was already cached.
This setting is optional for regeneration.

After the workflow completes successfully:

1. Copy the raw `.sscore` files from the PLINK2_SCORE task directories in `work/`.
2. Copy `aggregated_scores.txt.gz` from the SCORE_AGGREGATE task directory in `work/`.
3. Check that the outputs contain 24 chromosome files and 50 scored variant IDs.
4. Check that each chromosome has four samples, except Y, which has the two XY samples.
5. Check that the matching reports contain two unmatched MT rows.
6. Update the matching reports, logs, version report and provenance with any changed outputs.

Do not calculate a replacement reference with v3.
