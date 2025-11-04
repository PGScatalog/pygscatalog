import numpy as np
import pytest
from bgen import BgenWriter


@pytest.fixture(scope="session")
def big_bgen_n_samples():
    return 10_000


@pytest.fixture(scope="session")
def big_bgen_n_variants():
    return 500_000


@pytest.fixture(scope="session")
def really_big_bgen(tmp_path_factory, big_bgen_n_samples, big_bgen_n_variants):
    base_dir = tmp_path_factory.mktemp("big_file", numbered=False)
    bgen_path = base_dir / "test_large_file.bgen"
    bgi_path = base_dir / "test_large_file.bgen.bgi"

    # this file ~1GB
    with BgenWriter(bgen_path, n_samples=big_bgen_n_samples) as bfile:
        for i in range(big_bgen_n_variants):
            geno = np.random.rand(big_bgen_n_samples, 3).astype(np.float16)
            bfile.add_variant(
                varid=f"var{i}",
                rsid="rsid{i}",
                chrom=f"chr{i}",
                pos=i,
                alleles=["A", "G"],
                genotypes=geno,
                bit_depth=1,
            )

    sample_ids = [f"sample{i}" for i in range(big_bgen_n_samples)]

    return sample_ids, bgen_path, bgi_path
