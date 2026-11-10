"""LOAD rejects stale caches and preserves source archives when expanding them."""

import gzip
import json
import os
import shutil
import zipfile

import pytest
import zarr

from pgscatalog.calc.cli.load import (
    CACHE_VERSION,
    UNSUPPORTED_CACHE,
    source_metadata,
    unzip_zarr,
    validate_cached_source,
)

from .conftest import run_cli
from .test_sampleset_frequencies import _inputs


def _load(targets, scorefile, out, monkeypatch, cache=None, sample=None):
    argv = [
        "pgsc_calc",
        "load",
        "--target_genomes",
        *map(str, targets),
        "--sampleset",
        "test",
        "--format",
        "bgen" if sample else "vcf",
        "--score_paths",
        str(scorefile),
        "--cache_dir",
        str(out),
    ]
    if sample is not None:
        argv.extend(["--bgen_sample_file", str(sample)])
    if cache is not None:
        argv.extend(["--zarr_zip_file", str(cache)])
    run_cli(argv, monkeypatch)
    return out / "genotypes.zarr.zip"


def _attrs(archive, shard="test/ALL/whole.vcf.gz"):
    with zipfile.ZipFile(archive) as z:
        return json.loads(z.read(f"genotypes.zarr/{shard}/zarr.json"))["attributes"]


def test_cache_reuse_expands_without_mutating_source(tmp_path, monkeypatch):
    target, scorefile, header, _ = _inputs(tmp_path)
    first_score = tmp_path / "first.txt.gz"
    with gzip.open(scorefile, "rt") as source, gzip.open(first_score, "wt") as out:
        out.writelines(source.readlines()[:2])
    original = _load([target], first_score, tmp_path / "first", monkeypatch)
    before = original.read_bytes()
    attrs = _attrs(original)
    assert attrs["cache_version"] == CACHE_VERSION
    assert attrs["source_files"] == source_metadata(
        {"genotype": target, "index": target.with_name(target.name + ".csi")}
    )
    assert json.loads((tmp_path / "first/samples.json").read_text()) == {
        "shard": "test/ALL/whole.vcf.gz",
        "samples": list(header.samples),
    }

    # Moving unchanged inputs is allowed; preserving timestamps is essential.
    moved = tmp_path / "moved"
    moved.mkdir()
    for path in (target, target.with_name(target.name + ".csi")):
        shutil.copy2(path, moved / path.name)
    expanded = _load(
        [moved / target.name], scorefile, tmp_path / "expanded", monkeypatch, original
    )
    assert original.read_bytes() == before
    assert _attrs(expanded) == attrs
    unzip_zarr(expanded, tmp_path / "inspect")
    group = zarr.open_group(tmp_path / "inspect/genotypes.zarr/test/ALL/whole.vcf.gz")
    assert list(group["meta/chr_name"][...]) == ["1", "2"]
    assert list(group["meta/allele_counts"][...]) == ["20/80", "80/20"]
    reused = _load(
        [moved / target.name], scorefile, tmp_path / "reused", monkeypatch, expanded
    )
    assert _attrs(reused) == attrs


@pytest.mark.parametrize("role", ["genotype", "index"])
def test_changed_input_rejects_reuse(tmp_path, monkeypatch, role):
    target, scorefile, _, _ = _inputs(tmp_path)
    cache = _load([target], scorefile, tmp_path / "first", monkeypatch)
    changed = target if role == "genotype" else target.with_name(target.name + ".csi")
    stat = changed.stat()
    os.utime(changed, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1000))
    with pytest.raises(ValueError, match="Cache input metadata mismatch") as err:
        _load([target], scorefile, tmp_path / "second", monkeypatch, cache)
    assert "test/ALL/whole.vcf.gz" in str(err.value)
    assert str(changed) in str(err.value)
    assert "recorded=" in str(err.value) and "current=" in str(err.value)
    assert "publish a replacement" in str(err.value)


@pytest.mark.parametrize("change", ["version", "missing_metadata", "missing_leaf"])
def test_old_or_wrong_shard_cache_rejected(tmp_path, monkeypatch, change):
    target, scorefile, _, _ = _inputs(tmp_path)
    cache = _load([target], scorefile, tmp_path / "first", monkeypatch)
    incompatible = tmp_path / "old.zip"
    leaf = "genotypes.zarr/test/ALL/whole.vcf.gz/zarr.json"
    with zipfile.ZipFile(cache) as source, zipfile.ZipFile(incompatible, "w") as out:
        for name in source.namelist():
            data = source.read(name)
            if name == leaf:
                if change == "missing_leaf":
                    continue
                data = json.loads(data)
                if change == "version":
                    data["attributes"]["cache_version"] = 0
                else:
                    del data["attributes"]["source_files"]
                data = json.dumps(data).encode()
            out.writestr(name, data)
    with pytest.raises(ValueError) as err:
        _load([target], scorefile, tmp_path / "second", monkeypatch, incompatible)
    assert str(err.value) == UNSUPPORTED_CACHE


def test_empty_load_still_emits_sample_manifest(tmp_path, monkeypatch):
    target, scorefile, header, _ = _inputs(tmp_path)
    cache = _load([f"{target}:3"], scorefile, tmp_path / "empty", monkeypatch)
    assert not cache.exists()
    assert json.loads((tmp_path / "empty/samples.json").read_text()) == {
        "shard": "test/3/whole.vcf.gz",
        "samples": list(header.samples),
    }


def test_cannot_overwrite_supplied_cache(tmp_path, monkeypatch):
    target, scorefile, _, _ = _inputs(tmp_path)
    cache = _load([target], scorefile, tmp_path / "first", monkeypatch)
    before = cache.read_bytes()
    with pytest.raises(ValueError, match="preserve the supplied archive"):
        _load([target], scorefile, cache.parent, monkeypatch, cache)
    assert cache.read_bytes() == before


def test_optional_sample_file_metadata_is_validated(tmp_path):
    paths = {role: tmp_path / role for role in ("genotype", "index", "sample")}
    for path in paths.values():
        path.write_text("original")
    metadata = source_metadata(paths)
    cache = tmp_path / "cache.zip"
    with zipfile.ZipFile(cache, "w") as out:
        out.writestr(
            "genotypes.zarr/test/1/input.bgen/zarr.json",
            json.dumps(
                {
                    "attributes": {
                        "cache_version": CACHE_VERSION,
                        "source_files": metadata,
                    }
                }
            ),
        )
    validate_cached_source(cache, "test/1/input.bgen", paths, metadata)
    paths["sample"].write_text("changed size")
    with pytest.raises(ValueError, match="sample file.*recorded=.*current="):
        validate_cached_source(
            cache, "test/1/input.bgen", paths, source_metadata(paths)
        )


def test_multi_target_load_emits_each_shards_samples(tmp_path, monkeypatch):
    target, scorefile, header, _ = _inputs(tmp_path)
    _load([f"{target}:1", f"{target}:2"], scorefile, tmp_path / "multi", monkeypatch)
    assert json.loads((tmp_path / "multi/samples.json").read_text()) == {
        "shards": [
            {"shard": f"test/{chrom}/whole.vcf.gz", "samples": list(header.samples)}
            for chrom in ("1", "2")
        ]
    }


def test_bgen_load_records_and_validates_original_sample_and_index(
    tmp_path, monkeypatch, unphased_bgen_path, bgen_sample
):
    sample = tmp_path / bgen_sample.name
    shutil.copy2(bgen_sample, sample)
    scorefile = tmp_path / "score.txt.gz"
    with gzip.open(scorefile, "wt") as out:
        out.write(
            "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\t"
            "effect_type\tis_duplicated\taccession\trow_nr\n"
            "1\t44465567\tA\tG\t1\tadditive\tFalse\tPGS999998\t0\n"
        )
    target = f"{unphased_bgen_path}:1"
    paths = {
        "genotype": unphased_bgen_path,
        "index": unphased_bgen_path.with_name(unphased_bgen_path.name + ".bgi"),
        "sample": sample,
    }
    original = source_metadata(paths)
    cache = _load([target], scorefile, tmp_path / "bgen", monkeypatch, sample=sample)
    assert _attrs(cache, "test/1/tiny1000G.bgen")["source_files"] == original
    assert source_metadata(paths) == original
    _load([target], scorefile, tmp_path / "reuse", monkeypatch, cache, sample)
    sample.write_text("changed and no longer a valid sample file")
    with pytest.raises(ValueError, match="sample file.*recorded=.*current="):
        _load([target], scorefile, tmp_path / "changed", monkeypatch, cache, sample)
