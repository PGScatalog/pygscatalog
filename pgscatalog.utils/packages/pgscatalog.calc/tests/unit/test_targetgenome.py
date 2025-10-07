import pytest

from pgscatalog.calc import GenomeFileType, TargetGenome


@pytest.mark.parametrize(
    "filetype,path_fixture,sample_file_fixture",
    [
        (GenomeFileType.BGEN, "unphased_bgen_path", "bgen_sample"),
        (GenomeFileType.VCF, "vcf_path", None),
    ],
)
def test_targetgenome(
    filetype,
    path_fixture,
    sample_file_fixture,
    request,
    test_positions,
    tmp_path_factory,
    sampleset_name,
    sample_ids,
):
    cache_dir = tmp_path_factory.mktemp("cache")
    target_path = request.getfixturevalue(path_fixture)
    kwargs = {
        "target_path": target_path,
        "cache_dir": cache_dir,
        "sampleset": sampleset_name,
    }
    if sample_file_fixture:
        kwargs["sample_file"] = request.getfixturevalue(sample_file_fixture)

    target = TargetGenome(**kwargs)

    # test the public API of TargetGenome
    assert target.filetype == filetype
    assert target.variants_db_path.parent == cache_dir
    assert target.sampleset == sampleset_name
    assert target.cache_dir == cache_dir
    assert target.samples == sample_ids
    assert target.zarr_group.path == f"{sampleset_name}/{target_path.name}"
    assert target.chrom is None
    assert target.filename == target_path.name

    # cache logic and side effects are checked in test_cache.py
    assert target.cached_positions == frozenset()
    target.cache_variants(test_positions)
    assert target.cached_positions == frozenset(test_positions)
