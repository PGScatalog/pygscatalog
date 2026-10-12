"""Acceptance tests for the nested /{sampleset}/{chromosome}/{filename} cache
hierarchy (ADR 0003).

Scoring these shapes lives in test_four_stage_pipeline.py; this file is only about
where things land in the cache.
"""

import zipfile

import zarr
import zarr.storage

from pgscatalog.calc.cli.load import unzip_zarr

from .conftest import SAMPLESET, load_target, load_targets, run_cli


def _leaf_prefix(chrom: str, filename: str) -> str:
    return f"genotypes.zarr/{SAMPLESET}/{chrom}/{filename}/"


def test_two_files_per_chromosome(grid_targets, scorefile, tmp_path, monkeypatch):
    """Each shard is stored at its own chromosome, so basenames no longer collide"""
    zips = load_targets(grid_targets, tmp_path, scorefile, monkeypatch)

    for (path, chrom), zip_path in zip(grid_targets, zips, strict=True):
        names = zipfile.ZipFile(zip_path).namelist()
        prefix = _leaf_prefix(chrom, path.name)
        assert any(n.startswith(prefix + "genotypes/") for n in names), names[:5]


def test_unsplit_input_uses_all_sentinel(
    unsplit_target, scorefile, tmp_path, monkeypatch
):
    """Whole-genome input lands under the ALL sentinel chromosome"""
    zip_path = load_target(
        str(unsplit_target), tmp_path / "cache", scorefile, monkeypatch
    )

    names = zipfile.ZipFile(zip_path).namelist()
    prefix = _leaf_prefix("ALL", unsplit_target.name)
    assert any(n.startswith(prefix + "genotypes/") for n in names), names[:5]


def test_two_whole_genome_batches(
    sample_batch_targets, scorefile, tmp_path, monkeypatch
):
    """Two whole-genome VCFs holding disjoint samples are sample batches, not an
    ambiguous samplesheet: they load together and land side by side under ALL"""
    # one LOAD, both inputs: the old unsplit-chromosome rule rejected this
    run_cli(
        [
            "pgsc_calc",
            "load",
            "--target_genomes",
            *[str(p) for p in sample_batch_targets],
            "--sampleset",
            SAMPLESET,
            "--format",
            "vcf",
            "--cache_dir",
            str(tmp_path / "cache"),
            "--score_paths",
            str(scorefile),
        ],
        monkeypatch,
    )

    names = zipfile.ZipFile(tmp_path / "cache" / "genotypes.zarr.zip").namelist()
    for path in sample_batch_targets:
        prefix = _leaf_prefix("ALL", path.name) + "genotypes/"
        assert any(n.startswith(prefix) for n in names), names[:5]


def test_partial_extract_of_one_leaf(
    chromosome_targets, scorefile, tmp_path, monkeypatch
):
    """Extracting a single leaf brings every ancestor zarr.json with it, whatever
    the depth of the hierarchy"""
    path, chrom = chromosome_targets[0]
    zip_path = load_target(
        f"{path}:{chrom}", tmp_path / "cache", scorefile, monkeypatch
    )

    extract_dir = tmp_path / "extracted"
    extract_dir.mkdir()
    group_path = f"genotypes.zarr/{SAMPLESET}/{chrom}/{path.name}"
    unzip_zarr(zip_path=zip_path, cache_path=extract_dir, group_path=group_path)

    # the hierarchy only parses if root, sampleset and chromosome metadata is present
    store = zarr.storage.LocalStore(root=extract_dir / "genotypes.zarr")
    root = zarr.open_group(store=store, mode="r")
    leaf = root[f"{SAMPLESET}/{chrom}/{path.name}"]
    assert "genotypes" in leaf

    # and nothing else got dragged along
    assert list(root[SAMPLESET].group_keys()) == [chrom]


def test_absence_semantics_lands_on_the_leaf(
    unsplit_target, scorefile, tmp_path, monkeypatch
):
    """LOAD is fanned out per shard and merging zips is last-write-wins on an
    ancestor zarr.json, so the declaration goes on the leaf where a disagreement
    stays detectable"""
    cache_dir = tmp_path / "cache"
    load_target(str(unsplit_target), cache_dir, scorefile, monkeypatch)

    extract_dir = tmp_path / "extracted"
    extract_dir.mkdir()
    unzip_zarr(cache_dir / "genotypes.zarr.zip", extract_dir)

    store = zarr.storage.LocalStore(root=extract_dir / "genotypes.zarr")
    root = zarr.open_group(store=store, mode="r")
    leaf = root[f"{SAMPLESET}/ALL/{unsplit_target.name}"]

    assert leaf.attrs["absence_semantics"] == "unobserved"
