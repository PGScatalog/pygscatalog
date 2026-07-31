"""Acceptance tests for split-apply-combine scoring.

Scoring runs as four stages that fan out at different levels of the cache
hierarchy: MATCH (per chromosome) -> MATCH_COMBINE (once) -> SCORE (per shard) ->
SCORE_COMBINE (once). Every target genome shape here holds the same
(sample, variant) pairs, so all of them must produce the same scores whatever the
cache is split by.

Float addition isn't associative and one da.dot over every variant becomes N
partial dots plus a duckdb SUM, so scores are compared with pytest.approx.
"""

import gzip
import pathlib

import polars as pl
import pysam
import pytest
import zarr
import zarr.storage

from pgscatalog.calc.cli.load import unzip_zarr
from pgscatalog.calc.lib.score.score import ScoreStage
from pgscatalog.calc.lib.score.stage import scratch_dir

from .conftest import (
    SAMPLESET,
    _subset_samples,
    _write_vcf,
    load_targets,
    run_cli,
    run_pipeline,
)

KEY = ["sampleset", "accession", "sample_id"]


def _assert_matches_golden(scores: pl.DataFrame, golden: pl.DataFrame) -> None:
    scores, golden = scores.sort(KEY), golden.sort(KEY)

    assert scores[KEY].to_dicts() == golden[KEY].to_dicts()
    assert scores["n_matched"].to_list() == golden["n_matched"].to_list()
    assert scores["allele_count"].to_list() == golden["allele_count"].to_list()
    for column in ("score", "dosage_sum", "score_avg"):
        assert scores[column].to_list() == pytest.approx(
            golden[column].to_list(), abs=1e-6
        ), column


def test_unsplit_input(unsplit_target, scorefile, golden_scores, tmp_path, monkeypatch):
    """Whole-genome input: one chromosome group (ALL), one shard.

    One shard means one dot product in the same order the golden file was
    generated with, so this shape compares exactly rather than approximately -
    which is what pins autosomal scores while haploid support is added around
    them (ADR 0007).
    """
    scores = run_pipeline([unsplit_target], tmp_path, scorefile, monkeypatch)
    _assert_matches_golden(scores, golden_scores)

    scores, golden = scores.sort(KEY), golden_scores.sort(KEY)
    for column in ("score", "dosage_sum", "score_avg"):
        assert scores[column].to_list() == golden[column].to_list(), column


def test_chr_prefixed_input(
    chr_prefixed_target, scorefile, golden_scores, tmp_path, monkeypatch
):
    """A chr-prefixed VCF scores identically to an unprefixed one over the same data.

    The VCF reader normalised nothing, so `bcftools view -R 1` returned no regions
    and every variant was reported unmatched rather than failing (ADR 0008).
    """
    scores = run_pipeline([chr_prefixed_target], tmp_path, scorefile, monkeypatch)
    _assert_matches_golden(scores, golden_scores)


def test_variant_axis_only(
    chromosome_targets, scorefile, golden_scores, tmp_path, monkeypatch
):
    """One file per chromosome"""
    scores = run_pipeline(chromosome_targets, tmp_path, scorefile, monkeypatch)
    _assert_matches_golden(scores, golden_scores)


def test_sample_axis_only(
    sample_batch_targets, scorefile, golden_scores, tmp_path, monkeypatch
):
    """Two whole-genome files holding disjoint samples.

    This is the assertion that catches the per-leaf dedup bug, and nothing else in
    the suite does: keeping one match row per scoring variant across the whole run
    would leave one batch with no weight matrix row, and its samples would vanish
    from the output with nothing logged.
    """
    scores = run_pipeline(sample_batch_targets, tmp_path, scorefile, monkeypatch)
    _assert_matches_golden(scores, golden_scores)


def test_chromosome_by_sample_batch_grid(
    grid_targets, scorefile, golden_scores, tmp_path, monkeypatch
):
    """Both axes at once: 2 chromosomes x 2 sample batches"""
    scores = run_pipeline(grid_targets, tmp_path, scorefile, monkeypatch)
    _assert_matches_golden(scores, golden_scores)


def test_n_matched_is_genome_wide(
    grid_targets, scorefile, n_score_variants, tmp_path, monkeypatch
):
    """n_matched is a count of scoring variants, never a per-shard count summed
    four times"""
    scores = run_pipeline(grid_targets, tmp_path, scorefile, monkeypatch)

    assert scores["n_matched"].unique().to_list() == [n_score_variants]
    assert scores.shape[0] == 3202  # every sample scored exactly once
    assert scores["score"].n_unique() > 1


def test_two_accessions_over_disjoint_variants(
    chromosome_targets, target_records, scorefile, tmp_path, monkeypatch
):
    """Two accessions covering different chromosomes, so each selects a different
    slice of the weight matrix.

    weight_mat_row_nr is used directly as a row index into the dosage array, so
    this is what catches it drifting out of step with the array's group-major
    layout - the score column stays correct either way, but dosage_sum and
    allele_count do not.
    """
    _, by_chrom = target_records
    chrom_scorefiles = []
    for i, (chrom, records) in enumerate(by_chrom.items(), start=1):
        path = tmp_path / f"PGS99999{i}.txt.gz"
        accession = f"PGS99999{i}"
        with gzip.open(path, "wt") as f:
            f.write(
                "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\t"
                "effect_type\tis_duplicated\taccession\trow_nr\n"
            )
            for row_nr, rec in enumerate(records):
                f.write(
                    f"{chrom}\t{rec.pos}\t{rec.alts[0]}\t{rec.ref}\t0.1\t"
                    f"additive\tFalse\t{accession}\t{row_nr}\n"
                )
        chrom_scorefiles.append((accession, path, len(records)))

    scores = run_pipeline(
        chromosome_targets,
        tmp_path / "work",
        scorefile,
        monkeypatch,
        score_paths=[p for _, p, _ in chrom_scorefiles],
    )

    # each accession only ever sees its own chromosome's variants, so a sample
    # can't count more alleles than that chromosome holds
    for accession, _, n_variants in chrom_scorefiles:
        rows = scores.filter(pl.col("accession") == accession)
        assert rows["n_matched"].unique().to_list() == [n_variants]
        assert rows["allele_count"].max() <= n_variants * 2
        assert rows["dosage_sum"].max() <= n_variants * 2


def test_shared_samples_under_one_chromosome_are_rejected(
    variant_split_targets, scorefile, tmp_path, monkeypatch
):
    """Splitting variants within a chromosome leaves two shards sharing every
    sample, which would double-count dosage"""
    zips = load_targets(variant_split_targets, tmp_path, scorefile, monkeypatch)
    chrom = variant_split_targets[0][1]
    same_chrom = [
        z for z, (_, c) in zip(zips, variant_split_targets, strict=True) if c == chrom
    ]

    with pytest.raises(ValueError) as excinfo:
        run_cli(
            [
                "pgsc_calc",
                "match",
                "--zarr_zip_file",
                *[str(z) for z in same_chrom],
                "--score_paths",
                str(scorefile),
                "--out_dir",
                str(tmp_path / "match"),
            ],
            monkeypatch,
        )

    message = str(excinfo.value)
    # both offending shards are named, so the user knows what to look at
    for path, _ in variant_split_targets[:2]:
        assert path.name in message
    assert "share" in message


def test_ragged_shards_are_rejected_without_suggesting_homref(
    target_records, scorefile, tmp_path, monkeypatch
):
    """Ragged shards are ambiguous in cause: dropped monomorphic sites (where
    homozygous reference is true) and a re-imputed sample subset (where it is
    false) produce the same shape. `homref` is a claim about how a callset was
    generated, not a remedy for an error message.
    """
    header, by_chrom = target_records
    chrom, records = next(iter(by_chrom.items()))
    samples = list(header.samples)

    # rectangular batches would be fine; drop a variant from one of them
    whole = tmp_path / "whole.vcf.gz"
    _write_vcf(whole, header, records)

    ragged = []
    for i, half in enumerate((samples[::2], samples[1::2]), start=1):
        path = _subset_samples(whole, tmp_path / f"batch{i}.vcf.gz", half)
        if i == 2:
            trimmed = tmp_path / "batch2.trimmed.vcf.gz"
            with pysam.VariantFile(str(path)) as vcf:
                _write_vcf(trimmed, vcf.header.copy(), list(vcf.fetch())[1:])
            path = trimmed
        ragged.append((path, chrom))

    zips = load_targets(ragged, tmp_path / "load", scorefile, monkeypatch)

    with pytest.raises(ValueError) as excinfo:
        run_cli(
            [
                "pgsc_calc",
                "match",
                "--zarr_zip_file",
                *[str(z) for z in zips],
                "--score_paths",
                str(scorefile),
                "--out_dir",
                str(tmp_path / "match"),
            ],
            monkeypatch,
        )

    message = str(excinfo.value)
    assert "different variants" in message
    assert "homref" not in message.lower()


def test_accession_below_min_overlap_is_rejected_genome_wide(
    chromosome_targets, target_records, tmp_path, monkeypatch
):
    """The gate is genome-wide, so an accession can sit above the threshold on one
    chromosome and still be rejected"""
    _, by_chrom = target_records
    chrom_1 = by_chrom["1"]

    # every variant on chr1 matches; every variant on chr8 is a fake position that
    # can't match, dragging the genome-wide rate below min_overlap
    path = tmp_path / "PGS999997.txt.gz"
    with gzip.open(path, "wt") as f:
        f.write(
            "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\t"
            "effect_type\tis_duplicated\taccession\trow_nr\n"
        )
        row_nr = 0
        for rec in chrom_1:
            f.write(
                f"1\t{rec.pos}\t{rec.alts[0]}\t{rec.ref}\t0.1\tadditive\tFalse\t"
                f"PGS999997\t{row_nr}\n"
            )
            row_nr += 1
        # enough unmatchable variants to put the rate under 75%
        for i in range(len(chrom_1)):
            f.write(
                f"8\t{9_000_000 + i}\tA\tC\t0.1\tadditive\tFalse\tPGS999997\t{row_nr}\n"
            )
            row_nr += 1

    with pytest.raises(ValueError, match="No scores passed matching threshold"):
        run_pipeline(
            chromosome_targets,
            tmp_path / "work",
            path,
            monkeypatch,
            min_overlap="0.75",
        )


def test_sparse_scoring_file_does_not_trip_the_union_check(
    chromosome_targets, target_records, tmp_path, monkeypatch
):
    """A scoring file covering only some chromosomes leaves the others empty in the
    cache, which is not a sample mismatch"""
    _, by_chrom = target_records
    path = tmp_path / "PGS999996.txt.gz"
    with gzip.open(path, "wt") as f:
        f.write(
            "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\t"
            "effect_type\tis_duplicated\taccession\trow_nr\n"
        )
        for row_nr, rec in enumerate(by_chrom["1"]):
            f.write(
                f"1\t{rec.pos}\t{rec.alts[0]}\t{rec.ref}\t0.1\tadditive\tFalse\t"
                f"PGS999996\t{row_nr}\n"
            )

    # only chromosome 1 gets loaded, because chromosome 8 has nothing to cache
    scores = run_pipeline([chromosome_targets[0]], tmp_path / "work", path, monkeypatch)

    assert scores["accession"].unique().to_list() == ["PGS999996"]
    assert scores.shape[0] == 3202
    assert scores["n_matched"].unique().to_list() == [len(by_chrom["1"])]


def test_partial_scores_are_written_per_shard(
    grid_targets, scorefile, tmp_path, monkeypatch
):
    """Each SCORE task emits its own shard's contribution, covering only that
    shard's samples"""
    run_pipeline(grid_targets, tmp_path, scorefile, monkeypatch)

    partials = sorted((tmp_path / "score").glob("*/partial_scores.parquet"))
    assert len(partials) == len(grid_targets)

    shards = []
    for partial in partials:
        df = pl.read_parquet(partial)
        assert df["sampleset"].unique().to_list() == [SAMPLESET]
        # sample batches hold half the cohort each
        assert df["sample_id"].n_unique() == 1601
        # each partial is tagged with the one shard it came from, which is how
        # SCORE_COMBINE re-asserts uniqueness across the process boundary
        assert df["zarr_group"].n_unique() == 1
        shards.append(df["zarr_group"][0])

    assert len(set(shards)) == len(grid_targets)


def test_a_shard_reaching_combine_twice_is_rejected(
    grid_targets, scorefile, tmp_path, monkeypatch
):
    """Partial scores are summed per sample, so a shard that arrives twice doubles
    its contribution. Parquet carries no keys, so SCORE_COMBINE has to re-assert
    the constraint the score_table primary key can't express across files.
    """
    run_pipeline(grid_targets, tmp_path, scorefile, monkeypatch)
    partials = sorted((tmp_path / "score").glob("*/partial_scores.parquet"))

    with pytest.raises(ValueError, match="reached SCORE_COMBINE more than once"):
        run_cli(
            [
                "pgsc_calc",
                "combine",
                "--partial_scores",
                *[str(p) for p in [*partials, partials[0]]],
                "--out_dir",
                str(tmp_path / "combine_again"),
            ],
            monkeypatch,
        )


def test_dosage_array_is_sized_to_one_shard(
    grid_targets, scorefile, n_score_variants, tmp_path, monkeypatch
):
    """Peak memory per SCORE task scales with one shard rather than the sampleset,
    which is the point of the exercise.

    The mechanism is the dosage array's shape: it is allocated as (variants in this
    shard, samples in this shard), not (all variants, whole cohort).
    """
    run_pipeline(grid_targets, tmp_path, scorefile, monkeypatch)

    shapes = []
    for score_dir in sorted((tmp_path / "score").iterdir()):
        # the dosage array is scratch, not output
        store = zarr.storage.LocalStore(score_dir / "scratch")
        dosage = zarr.open_group(store=store, mode="r")[f"pgs/{SAMPLESET}"]["dosage"]
        shapes.append(dosage.shape)

    n_variants, n_samples = zip(*shapes, strict=True)
    # 2 chromosomes x 2 sample batches: each task sees half the cohort, and only
    # its own chromosome's variants
    assert set(n_samples) == {1601}
    assert sum(n_variants) == n_score_variants * 2  # each chromosome scored twice
    assert max(n_variants) < n_score_variants


def test_score_rejects_a_cache_holding_more_than_one_shard(
    chromosome_targets, scorefile, tmp_path, monkeypatch
):
    """weight_mat_row_nr indexes this task's dosage array, so SCORE staging two
    shards at once would silently misalign it"""
    zips = load_targets(chromosome_targets, tmp_path / "load", scorefile, monkeypatch)
    match_dir = tmp_path / "match"
    run_cli(
        [
            "pgsc_calc",
            "match",
            "--zarr_zip_file",
            str(zips[0]),
            "--score_paths",
            str(scorefile),
            "--out_dir",
            str(match_dir),
        ],
        monkeypatch,
    )
    combine_dir = tmp_path / "match_combine"
    run_cli(
        [
            "pgsc_calc",
            "match-combine",
            "--matches",
            str(match_dir / "matches.parquet"),
            "--samples",
            str(match_dir / "samples.parquet"),
            "--score_paths",
            str(scorefile),
            "--sampleset",
            SAMPLESET,
            # only one chromosome was matched, so the genome-wide rate is ~50%
            "--min_overlap",
            "0.4",
            "--out_dir",
            str(combine_dir),
        ],
        monkeypatch,
    )

    # stage both shards into one SCORE working directory by hand
    out_dir = tmp_path / "score"
    out_dir.mkdir()
    for zip_path in zips:
        unzip_zarr(zip_path, scratch_dir(out_dir))

    stage = ScoreStage(max_memory_gb=4, threads=1, out_dir=out_dir)
    stage.load_scores(scorefile_paths=[scorefile])
    stage.import_match_tables([combine_dir / "match_table.parquet"])

    with pytest.raises(ValueError, match="SCORE runs per shard"):
        stage.calculate_scores()


def test_score_rejects_the_ungated_match_table(
    unsplit_target, scorefile, tmp_path, monkeypatch
):
    """MATCH's matches.parquet and MATCH_COMBINE's match_table.parquet are both
    "the match table", and scoring the ungated one would score accessions that
    failed min_overlap without a crash, a warning or a NaN.

    Driven through the CLI because that is the wiring a real mistake goes through.
    """
    (zip_path,) = load_targets(
        [unsplit_target], tmp_path / "load", scorefile, monkeypatch
    )
    match_dir = tmp_path / "match"
    run_cli(
        [
            "pgsc_calc",
            "match",
            "--zarr_zip_file",
            str(zip_path),
            "--score_paths",
            str(scorefile),
            "--out_dir",
            str(match_dir),
        ],
        monkeypatch,
    )

    score_dir = tmp_path / "score"
    with pytest.raises(ValueError) as excinfo:
        run_cli(
            [
                "pgsc_calc",
                "score",
                "--zarr_zip_file",
                str(zip_path),
                "--score_paths",
                str(scorefile),
                "--match_table",
                str(match_dir / "matches.parquet"),
                "--out_dir",
                str(score_dir),
            ],
            monkeypatch,
        )

    # a bare BinderException reads like a corrupt parquet, so the message has to
    # name both artefacts and say what would otherwise have been scored
    message = str(excinfo.value)
    assert "matches.parquet" in message
    assert "match_table.parquet" in message
    assert "min_overlap" in message

    # it failed before doing any scoring work
    assert not (score_dir / "partial_scores.parquet").exists()


def test_match_rejects_more_than_one_chromosome(
    chromosome_targets, scorefile, tmp_path, monkeypatch
):
    """MATCH is a per-chromosome task, so staging two chromosomes is a wiring
    mistake rather than a bigger job"""
    zips = load_targets(chromosome_targets, tmp_path / "load", scorefile, monkeypatch)

    with pytest.raises(ValueError, match="MATCH runs per chromosome"):
        run_cli(
            [
                "pgsc_calc",
                "match",
                "--zarr_zip_file",
                *[str(z) for z in zips],
                "--score_paths",
                str(scorefile),
                "--out_dir",
                str(tmp_path / "match"),
            ],
            monkeypatch,
        )


def test_published_outputs_are_enumerable_without_the_scratch(
    grid_targets, scorefile, tmp_path, monkeypatch
):
    """A Nextflow process declares its outputs by name, so they have to be
    distinguishable from the things a task only needed while it ran.

    The score database, the unzipped cache and SCORE's dosage arrays all go under
    `scratch`, so the top level of a task's directory is exactly what it publishes.
    """
    run_pipeline(grid_targets, tmp_path, scorefile, monkeypatch)

    published = {
        tmp_path / "match" / "1": {"matches.parquet", "samples.parquet"},
        tmp_path / "match" / "8": {"matches.parquet", "samples.parquet"},
        tmp_path / "match_combine": {"match_table.parquet", "logs"},
        tmp_path / "combine": {"scores.txt.gz"},
        # 2 chromosomes x 2 sample batches
        **{tmp_path / "score" / str(i): {"partial_scores.parquet"} for i in range(4)},
    }

    for out_dir, artefacts in published.items():
        assert {path.name for path in out_dir.iterdir()} == artefacts | {"scratch"}

    # and the scratch is where the unpublished things actually went, so the
    # assertions above aren't passing because nothing was written
    scratch = {path.name for path in (tmp_path / "score" / "0" / "scratch").iterdir()}
    assert {"genotypes.zarr", "pgs.db", "pgs"} <= scratch


def test_match_logs_are_written(grid_targets, scorefile, tmp_path, monkeypatch):
    """MATCH_COMBINE owns the match logs now that MATCH doesn't compute rates"""
    run_pipeline(grid_targets, tmp_path, scorefile, monkeypatch)

    logs = tmp_path / "match_combine" / "logs"
    summary = pl.read_csv(logs / "summary.csv")

    assert summary["is_match_rate_ok"].unique().to_list() == [True]
    assert any(pathlib.Path(logs).rglob("*.csv.gz"))
