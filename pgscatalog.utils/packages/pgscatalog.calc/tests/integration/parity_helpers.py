"""Shared stage driver and readers for committed v2 references."""

import csv
import gzip
import hashlib
import json
from collections import defaultdict
from decimal import Decimal

import polars as pl

from .conftest import SAMPLESET, run_cli


def run_case(
    targets,
    score_paths,
    out_dir,
    monkeypatch,
    *,
    match_flags=(),
    min_overlap="0.01",
    batch_size=None,
):
    """Run the public stages. Retain intermediate outputs for assertions."""

    def stage(command, *args):
        run_cli(["pgsc_calc", command, *map(str, args)], monkeypatch)

    batch = ["--batch_size", batch_size] if batch_size is not None else []
    archives, matches, samples, partials = [], [], [], []
    chromosomes = defaultdict(list)
    for index, (path, chrom) in enumerate(targets):
        cache = out_dir / "load" / str(index)
        stage(
            "load",
            "--target_genomes",
            f"{path}:{chrom}",
            "--sampleset",
            SAMPLESET,
            "--format",
            "vcf",
            "--score_paths",
            *score_paths,
            "--cache_dir",
            cache,
            *batch,
        )
        archive = cache / "genotypes.zarr.zip"
        archives.append(archive)
        chromosomes[chrom].append(archive)
    # MATCH receives every sample batch of a chromosome together, as the public stage
    # graph requires. SCORE processes each shard separately.
    for index, chromosome_archives in enumerate(chromosomes.values()):
        out = out_dir / "match" / str(index)
        stage(
            "match",
            "--zarr_zip_file",
            *chromosome_archives,
            "--score_paths",
            *score_paths,
            "--out_dir",
            out,
            *match_flags,
            *batch,
        )
        matches.append(out / "matches.parquet")
        samples.append(out / "samples.parquet")
    combined = out_dir / "matched"
    stage(
        "match-combine",
        "--matches",
        *matches,
        "--samples",
        *samples,
        "--score_paths",
        *score_paths,
        "--sampleset",
        SAMPLESET,
        "--min_overlap",
        min_overlap,
        "--out_dir",
        combined,
    )
    for index, archive in enumerate(archives):
        out = out_dir / "score" / str(index)
        stage(
            "score",
            "--zarr_zip_file",
            archive,
            "--score_paths",
            *score_paths,
            "--match_table",
            combined / "match_table.parquet",
            "--out_dir",
            out,
            *batch,
        )
        partials.append(out / "partial_scores.parquet")
    stage("combine", "--partial_scores", *partials, "--out_dir", out_dir / "combined")
    return pl.read_csv(out_dir / "combined" / "scores.txt.gz")


def read_reference(reference_dir):
    """Sum effect-type shards. Keep DENOM separate from ALLELE_CT."""
    paths = sorted(reference_dir.glob("*.sscore"))
    assert paths, reference_dir
    result = defaultdict(lambda: defaultdict(Decimal))
    for path in paths:
        # The filename is sampleset_chrom_effecttype_index.sscore. The sampleset can
        # contain underscores.
        chrom = path.stem.rsplit("_", 3)[1]
        with path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert rows, path
        assert len(rows) == len({r["IID"] for r in rows}), path
        score_columns = [
            name
            for name in rows[0]
            if name.endswith("_SUM") and name != "NAMED_ALLELE_DOSAGE_SUM"
        ]
        assert score_columns, path
        for row in rows:
            for column in score_columns:
                key = (chrom, column.removesuffix("_SUM"), row["IID"])
                for field, source in (
                    ("score", column),
                    ("allele_count", "ALLELE_CT"),
                    ("denom", "DENOM"),
                    ("dosage_sum", "NAMED_ALLELE_DOSAGE_SUM"),
                ):
                    value = Decimal(row[source])
                    assert value.is_finite(), (path, key, field)
                    result[key][field] += value
    return dict(result)


def validate_inputs(case_dir, reference_dir=None):
    """Check saved outputs against their inputs and both scoring-file formats."""
    reference_dir = reference_dir or case_dir / "pgsc_calc_v2.3.0"
    provenance = json.loads((reference_dir / "provenance.json").read_text())
    assert provenance["commit"] == "72ee54f2119d2fe65ac671327342dfc412bd0a3b"
    assert {p.name for p in reference_dir.glob("*.sscore")} == set(
        provenance["sscore_files"]
    )
    for name, digest in provenance["input_sha256"].items():
        assert hashlib.sha256((case_dir / name).read_bytes()).hexdigest() == digest, (
            name
        )
    fields = (
        "chr_name",
        "chr_position",
        "effect_allele",
        "other_allele",
        "effect_weight",
    )
    for pair in provenance["score_pairs"]:
        path = case_dir / pair["v3"]
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as handle:
            v3 = list(csv.DictReader(handle, delimiter="\t"))
        with (case_dir / pair["v2"]).open() as handle:
            v2 = list(
                csv.DictReader(
                    (line for line in handle if not line.startswith("#")),
                    delimiter="\t",
                )
            )
        assert [[r[f] for f in fields] for r in v3] == [
            [r[f] for f in fields] for r in v2
        ], pair
        assert {r["accession"] for r in v3} == {pair["v2"].removesuffix(".txt")}
        for processed, original in zip(v3, v2, strict=True):
            effect_type = "additive"
            if original.get("is_dominant", "").lower() in ("true", "1"):
                effect_type = "dominant"
            if original.get("is_recessive", "").lower() in ("true", "1"):
                effect_type = "recessive"
            assert processed["effect_type"] == effect_type, pair
