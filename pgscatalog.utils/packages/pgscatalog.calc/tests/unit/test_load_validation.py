"""Samplesheet validation that runs before any work (ADR 0003)"""

import argparse

import pytest

from pgscatalog.calc.cli._utils import absence_semantics
from pgscatalog.calc.cli.load import (
    AnnotatedGenomePath,
    check_chromosomes,
    check_leaf_paths,
    parse_genome_paths,
)
from pgscatalog.calc.lib.cache.targetgenome import TargetGenome

SAMPLESET = "test"


def _targets(*paths: str) -> list[AnnotatedGenomePath]:
    return [AnnotatedGenomePath.from_string(x) for x in paths]


def test_colliding_basenames_are_rejected():
    """Directories are stripped, so these two inputs share a zarr leaf"""
    with pytest.raises(ValueError, match=r"1/batch1\.vcf\.gz") as excinfo:
        check_leaf_paths(
            _targets("/runA/batch1.vcf.gz:1", "/runB/batch1.vcf.gz:1"), SAMPLESET
        )

    # both offending inputs are named, so the user knows what to rename
    assert "/runA/batch1.vcf.gz" in str(excinfo.value)
    assert "/runB/batch1.vcf.gz" in str(excinfo.value)


def test_colliding_unsplit_basenames_are_rejected():
    """Undeclared chromosomes collide under the ALL sentinel too"""
    with pytest.raises(ValueError, match=r"ALL/wgs\.vcf\.gz"):
        check_leaf_paths(_targets("/runA/wgs.vcf.gz", "/runB/wgs.vcf.gz"), SAMPLESET)


def test_same_basename_on_different_chromosomes_is_fine():
    check_leaf_paths(
        _targets("/run/batch1.vcf.gz:1", "/run/batch1.vcf.gz:2"), SAMPLESET
    )


@pytest.mark.parametrize("chrom", [None, "1"])
def test_validated_leaf_path_is_the_written_leaf_path(
    chrom, vcf_path, vcf_index_path, tmp_path, sampleset_name
):
    """The validated string and the written group are one string, including the ALL
    sentinel. If they drifted, the guard would pass and a colliding leaf would
    silently overwrite the other (ADR 0003)"""
    target = TargetGenome(
        target_path=vcf_path,
        target_index_path=vcf_index_path,
        cache_dir=tmp_path,
        sampleset=sampleset_name,
        chrom=chrom,
    )
    validated = AnnotatedGenomePath(path=str(vcf_path), chrom=chrom)

    # zarr_group.path is where the shard is written, not a second copy of the rule
    assert validated.leaf_path(sampleset_name) == target.zarr_group.path


def test_validation_runs_before_any_work(tmp_path):
    """Colliding inputs fail without touching the filesystem, so the paths don't
    have to exist"""
    with pytest.raises(ValueError):
        parse_genome_paths(
            target_genomes=["/runA/batch1.vcf.gz:1", "/runB/batch1.vcf.gz:1"],
            cache_dir=tmp_path,
            sampleset="test",
            bgen_sample_file=None,
        )

    assert not (tmp_path / "genotypes.zarr").exists()


def test_two_unsplit_genomes_are_accepted():
    """Two whole-genome files are two sample batches under ALL, not an error"""
    targets = _targets("/run/batch1.vcf.gz", "/run/batch2.vcf.gz")

    check_leaf_paths(targets, SAMPLESET)
    check_chromosomes(targets)  # only needs .chrom

    assert [x.leaf_path(SAMPLESET) for x in targets] == [
        f"{SAMPLESET}/ALL/batch1.vcf.gz",
        f"{SAMPLESET}/ALL/batch2.vcf.gz",
    ]


def test_absence_semantics_default_is_accepted():
    assert absence_semantics("unobserved") == "unobserved"


def test_homref_is_rejected_not_ignored():
    """Silently ignoring the declaration would score a sparse-called callset as
    if every absent variant were unknown"""
    with pytest.raises(argparse.ArgumentTypeError, match="not implemented yet"):
        absence_semantics("homref")


def test_unknown_absence_semantics_is_rejected():
    with pytest.raises(argparse.ArgumentTypeError, match="not valid"):
        absence_semantics("missing")
