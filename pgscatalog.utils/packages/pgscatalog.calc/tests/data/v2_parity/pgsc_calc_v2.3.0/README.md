# pgsc_calc v2.3.0 reference

These are outputs of the actual Nextflow workflow at upstream tag `v2.3.0`,
commit `72ee54f2119d2fe65ac671327342dfc412bd0a3b`, run with the release's pinned
Docker images. The workflow source was not modified. `provenance.json` records
image digests, the command, original input hashes and chromosome scoring scripts;
`software_versions.yml` is the workflow's own version report.

## Coverage and configuration

The pinned matcher (`pgscatalog.match` 0.4.0) filters targets to chromosomes
1–22, X and Y, so this reference covers **24 chromosomes and 50 variants**.
MT is excluded from its samplesheet. `PGS999997.txt` contains all 52 original
weights, positions and allele pairs, converted to the custom scoring-file
format expected by v2. The matching summary and full matching log show the two
unmatched MT rows. `--min_overlap 0.95` admits 50/52 matches, but not 49/52.
The integration test independently requires exactly 50 matched variants.
It compares selected matches by accession and scoring-file row, requiring the
same chromosome, position and matched effect allele. It also resolves v3's
effect-allele indices against the original VCF alleles. Unselected v2 candidates,
duplicate flags and ALT annotations are not compared: v2's ID-only join can
duplicate multiallelic candidates and attach the wrong ALT annotation without
changing the effect allele passed to PLINK.

The configuration explicitly supplies the committed sample-sex metadata during
VCF import and uses `--xchr-model 1` during scoring, matching the dosage semantics
of v3 and the standalone PLINK reference. `--split-par b38` changes no variants:
all X positions are non-PAR. `--keep_multiallelic` retains multiallelic variants.
The uncompressed `.sscore` format was requested to make review straightforward.
These are configured v2.3.0 results, not a claim of equivalence to v2 defaults.

`--only_score` skips ancestry and report generation. Resource limits are reduced
for the tiny input; scoring runs with one thread. Nextflow 25.10.0 and the arm64
images recorded in `provenance.json` were used.

## Committed outputs

- `synthetic_{chrom}_additive_0.sscore`: 24 original PLINK outputs emitted by
  the v2 workflow, unchanged. The test aggregates `PGS999997_SUM`, `ALLELE_CT`
  and `NAMED_ALLELE_DOSAGE_SUM` by `IID`, then derives the score average.
- `aggregated_scores.txt.gz`: the original output of v2's SCORE_AGGREGATE stage.
  The test verifies its SUM, DENOM and AVG against the raw chromosome results.
- `synthetic_summary.csv` and `synthetic_log.csv.gz`: original matching reports.
- `plink.log`: concatenated original scoring logs; `nextflow.log`: workflow
  console output showing successful completion.
- `PGS999997.txt`, `samplesheet.csv` and `reference.config`: exact run inputs.
  The saved samplesheet and configuration retain the original absolute paths.

The v3 comparison uses the same committed score file and the same 24 VCFs.
The separate standalone-PLINK case continues to cover all 25 chromosomes,
including MT. Neither test launches PLINK, Nextflow or Docker.

## Regeneration

From `pgscatalog.utils`, with Nextflow and Docker available, prepare a fresh
run directory and adapt only the saved absolute fixture path:

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

Use a directory shared with the Docker VM; on this host `/private/tmp` was not
shared. The recorded run used `NXF_OFFLINE=true` because its Nextflow plugin
was already cached. That setting is optional for regeneration.

After successful completion, copy the raw `.sscore` files from the PLINK2_SCORE
task directories and `aggregated_scores.txt.gz` from SCORE_AGGREGATE in `work/`.
Verify 24 chromosome files, 50 scored variant IDs, four samples except for the
two XY samples on Y, and two unmatched MT score rows. Refresh the matching
reports, logs, version report and provenance along with any changed outputs.
Do not derive a replacement reference from v3.
