"""`pgsc_calc run` executes every stage in sequence, and refuses a partial one.

The scores a serial run produces are checked everywhere else in this suite, which
drives `pgsc_calc run` for every target genome shape (ADR 0013). What's here is the
tree it publishes and the two ways it's allowed to stop: a stage that failed, and a
shard whose partial score never arrived.
"""

import polars as pl
import pytest

from pgscatalog.calc.cli import run, score

from .conftest import run_serial


def test_every_stage_a_serial_run_executes_is_declared():
    """Nothing detects a fifth stage added to main.nf, so the list `run` iterates
    is written down and pinned here: dropping one is then a choice rather than an
    oversight (ADR 0013).
    """
    assert run.STAGES == ("load", "match", "match-combine", "score", "combine")


def test_the_published_tree_holds_a_zip_per_shard_and_scores_at_the_root(
    chromosome_targets, scorefile, tmp_path, monkeypatch
):
    """LOAD runs once per target into its own cache, because SCORE takes one zip
    holding one shard. The final scores are copied out of the stage tree.
    """
    run_serial(chromosome_targets, tmp_path, scorefile, monkeypatch)

    stages = tmp_path / "stages"
    assert sorted(path.name for path in (stages / "load").iterdir()) == ["0", "1"]
    assert all(
        (stages / "load" / str(i) / "genotypes.zarr.zip").exists() for i in range(2)
    )
    assert (tmp_path / "scores.txt.gz").exists()


def test_a_failed_stage_aborts_the_run(
    chromosome_targets, scorefile, tmp_path, monkeypatch
):
    """SCORE_COMBINE sums what it's handed, so carrying on past a failed SCORE
    would report a quietly low score for every sample instead of an error.
    """

    def explode(_args):
        raise ValueError("SCORE fell over")

    # every stage is dispatched through the parser, which reads this name when it
    # builds SCORE's subparser
    monkeypatch.setattr(score, "score_cli", explode)

    with pytest.raises(ValueError, match="SCORE fell over"):
        run_serial(chromosome_targets, tmp_path, scorefile, monkeypatch)

    # SCORE_COMBINE never ran, and the completed stages are still on disk so the
    # failing one can be rerun by hand against them
    assert not (tmp_path / "stages" / "combine").exists()
    assert not (tmp_path / "scores.txt.gz").exists()
    assert (tmp_path / "stages" / "match_combine" / "match_table.parquet").exists()


def test_a_missing_partial_score_is_refused_by_name(
    chromosome_targets, scorefile, tmp_path, monkeypatch
):
    """A stage that exits cleanly without publishing is the case SCORE_COMBINE
    can't detect: it doesn't know how many shards existed, and `run` does.
    """
    score_cli = score.score_cli

    def lose_the_last_shard(args):
        score_cli(args)
        if args.out_dir.name == "1":
            (args.out_dir / "partial_scores.parquet").unlink()

    monkeypatch.setattr(score, "score_cli", lose_the_last_shard)

    with pytest.raises(FileNotFoundError, match="partial_scores.parquet"):
        run_serial(chromosome_targets, tmp_path, scorefile, monkeypatch)

    assert not (tmp_path / "stages" / "combine").exists()
    assert not (tmp_path / "scores.txt.gz").exists()


def test_a_stage_missing_from_the_declared_list_is_refused(
    chromosome_targets, scorefile, tmp_path, monkeypatch
):
    """The declared list is checked against what `run` actually dispatched, so it
    can't quietly stop describing the run.
    """
    monkeypatch.setattr(run, "STAGES", ("load", "match", "match-combine", "score"))

    with pytest.raises(ValueError, match="must execute every stage"):
        run_serial(chromosome_targets, tmp_path, scorefile, monkeypatch)


def test_serial_run_rejects_colliding_leaf_paths_before_loading(
    scorefile, tmp_path, monkeypatch
):
    targets = [
        (tmp_path / "first/same.vcf.gz", "1"),
        (tmp_path / "second/same.vcf.gz", "1"),
    ]
    with pytest.raises(ValueError, match="both cache to"):
        run_serial(targets, tmp_path / "work", scorefile, monkeypatch)
    assert not (tmp_path / "work/stages/load").exists()


def test_serial_run_skips_chromosomes_without_scoring_positions(
    chromosome_targets, scorefile, tmp_path, monkeypatch
):
    chrom = chromosome_targets[0][1]
    rows = pl.read_csv(scorefile, separator="\t").filter(
        pl.col("chr_name").cast(pl.String) == chrom
    )
    filtered = tmp_path / "one-chromosome.txt"
    rows.write_csv(filtered, separator="\t")
    scores = run_serial(chromosome_targets, tmp_path / "work", filtered, monkeypatch)
    assert scores["n_matched"].unique().to_list() == [rows.height]
    assert not (tmp_path / "work/stages/load/1").exists()
    with pytest.raises(ValueError, match="No target chromosomes"):
        run_serial(chromosome_targets[1:], tmp_path / "unused", filtered, monkeypatch)
