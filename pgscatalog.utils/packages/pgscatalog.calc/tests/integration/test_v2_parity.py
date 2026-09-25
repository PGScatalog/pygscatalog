"""All-chromosome score parity with the saved independent PLINK 2 reference.

PLINK is only used to regenerate the reference, never during this test.
See tests/data/v2_parity/README.md for the fixture and regeneration commands.
"""

from __future__ import annotations

import csv
import gzip
import pathlib
from collections import defaultdict
from decimal import Decimal

from .conftest import SAMPLESET, load_targets, run_cli

DATA_DIR = pathlib.Path(__file__).parent.parent / "data" / "v2_parity"
CHROMS = (*map(str, range(1, 23)), "X", "Y", "MT")
SCORE_FILE = DATA_DIR / "score.txt.gz"
REFERENCE = DATA_DIR / "plink_reference"


def test_all_chromosomes_match_plink_scores(tmp_path, monkeypatch):
    """Exercise LOAD through COMBINE with multiallelic and mixed-ploidy calls."""
    sample_ids = {"XX_1", "XX_2", "XY_1", "XY_2"}
    expected = defaultdict(Decimal)
    assert {p.name for p in REFERENCE.glob("*.sscore")} == {
        f"chr{chrom}.sscore" for chrom in CHROMS
    }
    for chrom in CHROMS:
        with (REFERENCE / f"chr{chrom}.sscore").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        wanted = {"XY_1", "XY_2"} if chrom == "Y" else sample_ids
        assert len(rows) == len(wanted), chrom
        assert {r["#IID"] for r in rows} == wanted, chrom
        for row in rows:
            contribution = Decimal(row["WEIGHT_SUM"])
            assert contribution.is_finite(), (chrom, row["#IID"])
            expected["PGS999997", row["#IID"]] += contribution
    # Round once after adding chromosome sums, never average WEIGHT_AVG.
    expected = {
        key: value.quantize(Decimal("0.000001")) for key, value in expected.items()
    }

    targets = [(DATA_DIR / f"target.chr{chrom}.vcf.gz", chrom) for chrom in CHROMS]
    zips = load_targets(
        targets,
        tmp_path / "load",
        SCORE_FILE,
        monkeypatch,
        sex_map=DATA_DIR / "sex.tsv",
    )

    def stage(command, *args):
        run_cli(["pgsc_calc", command, *map(str, args)], monkeypatch)

    # run does not expose --keep_multiallelic yet; exercise the public stages.
    matches, samples, partials = [], [], []
    for chrom, archive in zip(CHROMS, zips, strict=True):
        out = tmp_path / "match" / chrom
        stage(
            "match",
            "--zarr_zip_file",
            archive,
            "--score_paths",
            SCORE_FILE,
            "--keep_multiallelic",
            "--out_dir",
            out,
        )
        matches.append(out / "matches.parquet")
        samples.append(out / "samples.parquet")
    combined = tmp_path / "matched"
    stage(
        "match-combine",
        "--matches",
        *matches,
        "--samples",
        *samples,
        "--score_paths",
        SCORE_FILE,
        "--sampleset",
        SAMPLESET,
        # Current gate uses >. 51/52 < 0.99, so only all 52 can pass.
        "--min_overlap",
        "0.99",
        "--out_dir",
        combined,
    )
    for chrom, archive in zip(CHROMS, zips, strict=True):
        out = tmp_path / "score" / chrom
        stage(
            "score",
            "--zarr_zip_file",
            archive,
            "--score_paths",
            SCORE_FILE,
            "--match_table",
            combined / "match_table.parquet",
            "--out_dir",
            out,
        )
        partials.append(out / "partial_scores.parquet")
    stage("combine", "--partial_scores", *partials, "--out_dir", tmp_path / "combined")

    with gzip.open(tmp_path / "combined" / "scores.txt.gz", "rt") as handle:
        actual = list(csv.DictReader(handle))
    keys = {("PGS999997", sample) for sample in sample_ids}
    assert len(actual) == len(expected) == 4
    assert set(expected) == keys
    assert {(r["accession"], r["sample_id"]) for r in actual} == keys
    assert {r["sampleset"] for r in actual} == {SAMPLESET}
    assert all(int(r["n_matched"]) == 52 for r in actual)
    # Compare the public six-decimal score sums exactly, with no float tolerance.
    assert {
        (r["accession"], r["sample_id"]): Decimal(r["score"]) for r in actual
    } == expected
