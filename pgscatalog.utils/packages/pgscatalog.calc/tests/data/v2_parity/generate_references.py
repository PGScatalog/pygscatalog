"""Regenerate extended references with an unchanged v2.3.0 checkout and Docker.

Run deliberately outside pytest. The work directory must be new and shared with
Docker. This writes the selected case's reference outputs, never VCF inputs.
"""

import argparse
import csv
import hashlib
import json
import os
import pathlib
import shutil
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("case", choices=["arithmetic", "missingness", "matching", "zero_weight"])
p.add_argument(
    "--mode",
    choices=["default", "no_mean_imputation", "keep_ambiguous", "ignore_strand_flips"],
    default="default",
)
p.add_argument("--pipeline", type=pathlib.Path, required=True)
p.add_argument("--work-dir", type=pathlib.Path, required=True)
a = p.parse_args()
fixture = pathlib.Path(__file__).resolve().parent
case = fixture / a.case
run = a.work_dir.resolve()
run.mkdir(parents=True, exist_ok=False)
pipeline = a.pipeline.resolve()
commit = subprocess.check_output(
    ["git", "-C", str(pipeline), "rev-parse", "HEAD"], text=True
).strip()
if commit != "72ee54f2119d2fe65ac671327342dfc412bd0a3b":
    raise ValueError("Expected upstream v2.3.0 commit")
if subprocess.check_output(
    ["git", "-C", str(pipeline), "status", "--porcelain"], text=True
).strip():
    raise ValueError("Pipeline checkout must be clean")
images = json.loads((fixture / "pgsc_calc_v2.3.0/provenance.json").read_text())[
    "images"
]
for tag, digest in images.items():
    metadata = json.loads(subprocess.check_output(["docker", "image", "inspect", tag]))
    if not any(value.endswith("@" + digest) for value in metadata[0]["RepoDigests"]):
        raise ValueError(f"Image digest differs from the recorded release: {tag}")
with (run / "samplesheet.csv").open("w") as f:
    w = csv.writer(f)
    w.writerow(["sampleset", "path_prefix", "chrom", "format"])
    for vcf in sorted(case.glob("target.chr*.vcf.gz")):
        chrom = vcf.name.removeprefix("target.chr").removesuffix(".vcf.gz")
        w.writerow(["synthetic", str(vcf).removesuffix(".vcf.gz"), chrom, "vcf"])
args2 = "no-mean-imputation" if a.mode == "no_mean_imputation" else ""
(run / "reference.config").write_text(
    "params.custom_config_base = null\nparams.max_memory = '2.GB'\n"
    "params.max_cpus = 1\nexecutor.queueSize = 4\nprocess {\n"
    " withName: PLINK2_SCORE {\n"
    " ext.args = '--allow-extra-chr --xchr-model 1'\n ext.args2 = '"
    + args2
    + "'\n }\n}\n"
)
if a.mode == "ignore_strand_flips":
    with (run / "reference.config").open("a") as f:
        f.write(
            "process { withName: 'MATCH_VARIANTS|MATCH_COMBINE' { "
            "ext.args = '--ignore_strand_flips' } }\n"
        )
cmd = [
    "nextflow",
    "run",
    str(pipeline),
    "-profile",
    "docker",
    "-c",
    "reference.config",
    "--input",
    "samplesheet.csv",
    "--scorefile",
    str(case / "PGS*.txt"),
    "--target_build",
    "GRCh38",
    "--keep_multiallelic",
    "--min_overlap",
    "0.01",
    "--only_score",
    "--outdir",
    "results",
    "-ansi-log",
    "false",
]
if a.mode == "keep_ambiguous":
    cmd += ["--" + a.mode]
env = os.environ.copy()
with (run / "run.log").open("w") as f:
    result = subprocess.run(cmd, cwd=run, env=env, stdout=f, stderr=subprocess.STDOUT)
print("RUN", run, "EXIT", result.returncode, flush=True)
if result.returncode:
    print((run / "run.log").read_text()[-8000:])
    raise SystemExit(result.returncode)
ref = case / "pgsc_calc_v2.3.0"
if a.case == "matching" or a.mode != "default":
    ref /= a.mode
ref.mkdir(parents=True, exist_ok=True)
files = {}
for pattern in (
    "*.sscore",
    "aggregated_scores.txt.gz",
    "synthetic_log.csv.gz",
    "synthetic_summary.csv",
    "software_versions.yml",
):
    for path in (run / "work").rglob(pattern):
        if path.is_symlink():
            continue
        if path.name in files and files[path.name].read_bytes() != path.read_bytes():
            raise ValueError(("conflicting", path))
        files[path.name] = path
for name, path in files.items():
    shutil.copyfile(path, ref / name)
for name in ("samplesheet.csv", "reference.config", "run.log"):
    shutil.copyfile(run / name, ref / name)
inputs = [f for f in case.iterdir() if f.is_file() and f.suffix not in (".py", ".md")]
pairs = [
    {
        "v2": f.name,
        "v3": f.name + ".gz" if (case / (f.name + ".gz")).exists() else "score.txt.gz",
    }
    for f in sorted(case.glob("PGS*.txt"))
]
prov = {
    "release": "v2.3.0",
    "commit": "72ee54f2119d2fe65ac671327342dfc412bd0a3b",
    "command": cmd,
    "working_directory": str(run),
    "input_sha256": {
        f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs
    },
    "score_pairs": pairs,
    "sscore_files": sorted(name for name in files if name.endswith(".sscore")),
    "task_commands": {
        str(f.parent.relative_to(run / "work")): f.read_text()
        for f in (run / "work").rglob(".command.sh")
        if "plink2 " in f.read_text()
    },
    "images": images,
}
(ref / "provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
print("SAVED", ref, "FILES", list(files), flush=True)
