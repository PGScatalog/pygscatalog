"""The cache hierarchy: what a staged cache holds (ADR 0003)"""

import numpy as np
import pytest
import zarr
import zarr.storage

from pgscatalog.calc.lib.cache.genotypecache import GenotypeCache, meta_group_path

SAMPLESET = "test"


def _build_cache(path, leaves, sampleset=SAMPLESET):
    """A cache holding the given leaf groups, each carrying sample IDs

    Leaves are (path below the sampleset, samples) pairs.
    """
    store = zarr.storage.LocalStore(path)
    for leaf_path, samples in leaves:
        group = zarr.open_group(
            store=store, path=f"{sampleset}/{leaf_path}", mode="a", zarr_format=3
        )
        group.attrs["samples"] = list(samples)
        # a real leaf also holds meta and genotypes, which must not be mistaken
        # for shards: only leaves carry sample IDs
        zarr.open_group(
            store=store,
            path=meta_group_path(f"{sampleset}/{leaf_path}"),
            mode="a",
            zarr_format=3,
        )
    return GenotypeCache(path=path)


def _write_meta(path, leaf_path, positions, sampleset=SAMPLESET):
    """Give a leaf the 1D variant metadata LOAD writes"""
    meta = zarr.open_group(
        store=zarr.storage.LocalStore(path),
        path=meta_group_path(f"{sampleset}/{leaf_path}"),
        mode="a",
        zarr_format=3,
    )
    meta.create_array("chr_pos", data=np.array(positions, dtype="int64"))


def test_a_shards_chromosome_is_the_group_it_sits_under(tmp_path):
    """The chromosome is a shard's position in the tree, not a slice of its path
    (ADR 0003)"""
    cache = _build_cache(
        tmp_path / "genotypes.zarr",
        [("1/a.vcf.gz", ["s1", "s2"]), ("8/a.vcf.gz", ["s1", "s2"])],
    )

    assert {(s.path, s.chromosome) for s in cache.shards} == {
        ("test/1/a.vcf.gz", "1"),
        ("test/8/a.vcf.gz", "8"),
    }


def test_shards_are_found_at_whatever_depth_they_sit(tmp_path):
    """Leaves are the only nodes carrying sample IDs, so they're found by that
    rather than by counting levels"""
    cache = _build_cache(
        tmp_path / "genotypes.zarr", [("1/nested/deeper/a.vcf.gz", ["s1"])]
    )

    (shard,) = cache.shards
    assert shard.path == "test/1/nested/deeper/a.vcf.gz"
    assert shard.chromosome == "1"


def test_unsplit_input_sits_under_the_all_sentinel(tmp_path):
    cache = _build_cache(tmp_path / "genotypes.zarr", [("ALL/whole.vcf.gz", ["s1"])])

    assert cache.chromosome == "ALL"


def test_samples_and_absence_semantics_are_declared_per_shard(tmp_path):
    """Absence semantics defaults to unobserved when LOAD declared nothing
    (ADR 0004)"""
    cache = _build_cache(tmp_path / "genotypes.zarr", [("1/a.vcf.gz", ["s1", "s2"])])

    assert cache.shard.samples == ["s1", "s2"]
    assert cache.shard.absence_semantics == "unobserved"


def test_sample_ids_come_from_the_one_staged_shard(tmp_path):
    """What SCORE asks for: one leaf, so nothing to disagree with (ADR 0003)"""
    cache = _build_cache(tmp_path / "genotypes.zarr", [("1/a.vcf.gz", ["s1", "s2"])])

    assert cache.sample_ids == ["s1", "s2"]


def test_sample_ids_reject_shards_that_disagree(tmp_path):
    """Kept exactly as it was when one process saw every shard at once: the
    assumption stopped being load-bearing rather than being fixed (ADR 0003)"""
    cache = _build_cache(
        tmp_path / "genotypes.zarr",
        [("1/a.vcf.gz", ["s1", "s2"]), ("1/b.vcf.gz", ["s3", "s4"])],
    )

    with pytest.raises(
        NotImplementedError, match="Different samples across target genome files"
    ):
        _ = cache.sample_ids


def test_sample_ids_reject_a_cache_with_no_shards(tmp_path):
    """A traversal that finds no shards is an error, not an empty cache"""
    path = tmp_path / "genotypes.zarr"
    zarr.open_group(
        store=zarr.storage.LocalStore(path),
        path=f"{SAMPLESET}/1",
        mode="a",
        zarr_format=3,
    )

    with pytest.raises(TypeError, match="Couldn't get samples from zarr attribute"):
        _ = GenotypeCache(path=path).sample_ids


def test_variant_metadata_is_read_one_shard_at_a_time(tmp_path):
    """MATCH holds every leaf's `meta` for a chromosome, keyed by shard path"""
    path = tmp_path / "genotypes.zarr"
    cache = _build_cache(
        path, [("1/a.vcf.gz", ["s1", "s2"]), ("1/b.vcf.gz", ["s3", "s4"])]
    )
    _write_meta(path, "1/a.vcf.gz", [1, 2])
    _write_meta(path, "1/b.vcf.gz", [3])

    assert {
        shard_path: list(meta["chr_pos"])
        for shard_path, meta in cache.variant_metadata()
    } == {"test/1/a.vcf.gz": [1, 2], "test/1/b.vcf.gz": [3]}


def test_a_cache_holding_one_shard_is_what_score_wants(tmp_path):
    cache = _build_cache(tmp_path / "genotypes.zarr", [("1/a.vcf.gz", ["s1"])])

    assert cache.shard.path == "test/1/a.vcf.gz"


def test_score_rejects_a_cache_holding_more_than_one_shard(tmp_path):
    """weight_mat_row_nr indexes this task's dosage array, so staging two shards
    would silently misalign it"""
    cache = _build_cache(
        tmp_path / "genotypes.zarr",
        [("1/a.vcf.gz", ["s1"]), ("1/b.vcf.gz", ["s2"])],
    )

    with pytest.raises(ValueError) as excinfo:
        _ = cache.shard

    message = str(excinfo.value)
    assert "SCORE runs per shard" in message
    assert "test/1/a.vcf.gz" in message and "test/1/b.vcf.gz" in message


def test_match_rejects_a_cache_holding_more_than_one_chromosome(tmp_path):
    cache = _build_cache(
        tmp_path / "genotypes.zarr",
        [("1/a.vcf.gz", ["s1"]), ("8/a.vcf.gz", ["s1"])],
    )

    with pytest.raises(ValueError, match="MATCH runs per chromosome"):
        _ = cache.chromosome


def test_the_sampleset_is_read_from_the_store(tmp_path):
    cache = _build_cache(tmp_path / "genotypes.zarr", [("1/a.vcf.gz", ["s1"])])

    assert cache.sampleset == "test"


def test_more_than_one_sampleset_is_rejected(tmp_path):
    """Shards of one sampleset are slices of one callset, so two samplesets in one
    store means two cohorts were staged together"""
    path = tmp_path / "genotypes.zarr"
    _build_cache(path, [("1/a.vcf.gz", ["s1"])], sampleset="one")
    cache = _build_cache(path, [("1/a.vcf.gz", ["s2"])], sampleset="two")

    with pytest.raises(NotImplementedError, match="One than one sampleset"):
        _ = cache.sampleset


def test_a_store_with_no_sampleset_is_rejected(tmp_path):
    path = tmp_path / "genotypes.zarr"
    zarr.open_group(store=zarr.storage.LocalStore(path), mode="a", zarr_format=3)

    with pytest.raises(FileNotFoundError, match="can't find sampleset group"):
        _ = GenotypeCache(path=path).sampleset


def test_a_declared_sampleset_needs_no_store(tmp_path):
    """MATCH_COMBINE and SCORE_COMBINE have no genotypes.zarr to read it from"""
    cache = GenotypeCache(path=tmp_path / "genotypes.zarr", sampleset="test")

    assert cache.sampleset == "test"
