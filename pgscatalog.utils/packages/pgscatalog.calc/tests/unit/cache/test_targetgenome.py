import pytest
import os

from pgscatalog.calc.lib.cache.genomefiletypes import GenomeFileType
from pgscatalog.calc.lib.cache.targetgenome import TargetGenome


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

    # __repr__
    repr_str = repr(target)
    assert "TargetGenome" in repr_str
    assert f"target_path={os.fspath(target_path)!r}" in repr_str
    assert "chrom=None" in repr_str
    assert f"cache_dir={os.fspath(cache_dir)!r}" in repr_str

    # test the properties of TargetGenome
    assert target.filetype == filetype
    assert target.cache_dir == cache_dir
    assert target.sampleset == sampleset_name
    assert target.samples == sample_ids
    assert target.target_path == target_path
    assert target.filename == target_path.name
    assert target.chrom is None

    # zarr layout
    assert target.zarr_group.path == f"{sampleset_name}/{target_path.name}"
    assert target._zarr_archive_name == cache_dir / "genotypes.zarr"
    assert target._zarr_archive_name.exists()

    # caching variants
    target.cache_variants(test_positions)
    cached = target.cached_positions
    assert {"chr_name", "chr_pos"}.issubset(set(cached.columns))
    target.cache_variants(test_positions)
