"""Samplesheet validation that runs before any work (docs/tickets/0002, ADR 0003)"""

import pytest

from pgscatalog.calc.cli.load import (
    AnnotatedGenomePath,
    check_chromosomes,
    check_leaf_paths,
    parse_genome_paths,
)


def _targets(*paths: str) -> list[AnnotatedGenomePath]:
    return [AnnotatedGenomePath.from_string(x) for x in paths]


def test_colliding_basenames_are_rejected():
    """Directories are stripped, so these two inputs share a zarr leaf"""
    with pytest.raises(ValueError, match=r"1/batch1\.vcf\.gz") as excinfo:
        check_leaf_paths(_targets("/runA/batch1.vcf.gz:1", "/runB/batch1.vcf.gz:1"))

    # both offending inputs are named, so the user knows what to rename
    assert "/runA/batch1.vcf.gz" in str(excinfo.value)
    assert "/runB/batch1.vcf.gz" in str(excinfo.value)


def test_colliding_unsplit_basenames_are_rejected():
    """Undeclared chromosomes collide under the ALL sentinel too"""
    with pytest.raises(ValueError, match=r"ALL/wgs\.vcf\.gz"):
        check_leaf_paths(_targets("/runA/wgs.vcf.gz", "/runB/wgs.vcf.gz"))


def test_same_basename_on_different_chromosomes_is_fine():
    check_leaf_paths(_targets("/run/batch1.vcf.gz:1", "/run/batch1.vcf.gz:2"))


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

    check_leaf_paths(targets)
    check_chromosomes(targets)  # only needs .chrom

    assert [x.leaf_path for x in targets] == ["ALL/batch1.vcf.gz", "ALL/batch2.vcf.gz"]
