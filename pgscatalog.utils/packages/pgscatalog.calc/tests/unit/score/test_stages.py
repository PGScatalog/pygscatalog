"""Each stage of split-apply-combine takes only what it needs.

`ScorePipeline` was one module whose interface was the union of four stages' needs,
with three parts that were conditional on the stage and invisible in a type
signature. These are the assertions that the union is gone.
"""

import importlib

import pytest

from pgscatalog.calc.lib.score.match import MatchStage
from pgscatalog.calc.lib.score.match_combine import MatchCombineStage
from pgscatalog.calc.lib.score.score import ScoreStage
from pgscatalog.calc.lib.score.score_combine import ScoreCombineStage

# every public member of a stage is one its own CLI calls on every run, so there is
# nothing left that raises depending on which stage constructed it
STAGE_MEMBERS = {
    MatchStage: {
        "load_scores",
        "match_variants",
        "export_match_table",
        "export_samples",
    },
    MatchCombineStage: {
        "load_scores",
        "import_match_tables",
        "combine_matches",
        "export_full_match_log",
        "export_summary_match_log",
    },
    ScoreStage: {
        "load_scores",
        "import_match_tables",
        "calculate_scores",
        "export_partial_scores",
    },
    ScoreCombineStage: {"combine_scores", "export_scores"},
}


@pytest.mark.parametrize(("stage", "members"), STAGE_MEMBERS.items())
def test_a_stage_has_no_member_it_cannot_use(stage, members):
    assert {name for name in vars(stage) if not name.startswith("_")} == members


@pytest.mark.parametrize(
    ("module_name", "stage"),
    [
        ("pgscatalog.calc.cli.match", MatchStage),
        ("pgscatalog.calc.cli.match_combine", MatchCombineStage),
        ("pgscatalog.calc.cli.score", ScoreStage),
        ("pgscatalog.calc.cli.combine", ScoreCombineStage),
    ],
)
def test_each_cli_constructs_only_its_own_stage(module_name, stage):
    """The stage graph stays in Nextflow (ADR 0001), so a subcommand is one stage and
    no entrypoint can reach another stage's members"""
    module = importlib.import_module(module_name)
    constructed = {
        obj
        for obj in vars(module).values()
        if isinstance(obj, type) and obj in STAGE_MEMBERS
    }
    assert constructed == {stage}


def test_match_combine_cannot_be_constructed_without_a_sampleset(tmp_path):
    """It has no cache to read the sampleset from, so the stage that has to be told
    is the one that asks for it"""
    with pytest.raises(TypeError):
        MatchCombineStage(out_dir=tmp_path, max_memory_gb=1, threads=1)


def test_only_score_takes_an_impute_floor(tmp_path):
    """Nothing else calculates dosage, so nothing else has anything to say about it"""
    stage = ScoreStage(
        out_dir=tmp_path, max_memory_gb=1, threads=1, minimum_samples_for_impute=3
    )
    assert stage._n_minimum_samples_for_impute == 3

    for other in (MatchStage, ScoreCombineStage):
        with pytest.raises(TypeError):
            other(
                out_dir=tmp_path,
                max_memory_gb=1,
                threads=1,
                minimum_samples_for_impute=3,
            )
