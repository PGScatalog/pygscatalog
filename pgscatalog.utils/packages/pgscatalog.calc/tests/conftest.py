import csv
import gzip
import pathlib
import random

import dask.array as da
import numpy as np
import pysam
import pysam.bcftools
import pytest
import zarr

from pgscatalog.calc import Scorefiles, TargetGenome
from pgscatalog.calc.lib._dosage import calculate_effect_allele_dosage
from pgscatalog.calc.lib.constants import MISSING_GENOTYPE_SENTINEL_VALUE


@pytest.fixture
def test_positions() -> list[tuple[str, int]]:
    """Some positions in tiny1000G"""
    return [("8", 116938529), ("1", 44465567), ("10", 97611390)]


@pytest.fixture
def unphased_bgen_path() -> pathlib.Path:
    return (
        pathlib.Path(__file__).parent / "data" / "bgen" / "unphased" / "tiny1000G.bgen"
    )


@pytest.fixture
def phased_bgen_path() -> pathlib.Path:
    return pathlib.Path(__file__).parent / "data" / "bgen" / "phased" / "tiny1000G.bgen"


@pytest.fixture
def bgen_sample() -> pathlib.Path:
    return (
        pathlib.Path(__file__).parent / "data" / "bgen" / "phased" / "tiny1000G.sample"
    )

@pytest.fixture
def bgen_samples(bgen_sample) -> list[str]:
    with open(bgen_sample, "rt") as f:
        csvreader = csv.DictReader(f, delimiter=" ")
        return [row["ID_2"] for row in csvreader][1:]

@pytest.fixture(scope="session")
def vcf_path() -> pathlib.Path:
    """Path to an indexed version of 1000 Genomes, containing ~1000 variants randomly
    sampled from PGS001229"""
    return pathlib.Path(__file__).parent / "data" / "vcf" / "tiny1000G.vcf.gz"


@pytest.fixture(scope="session")
def pgs001229() -> pathlib.Path:
    return (
        pathlib.Path(__file__).parent
        / "data"
        / "structured"
        / "PGS001229_hmPOS_GRCh38.txt.gz"
    )


@pytest.fixture(scope="session")
def sampleset_name() -> str:
    return "test"


@pytest.fixture(scope="session")
def pgs000001() -> pathlib.Path:
    return (
        pathlib.Path(__file__).parent
        / "data"
        / "structured"
        / "PGS000001_hmPOS_GRCh38.txt.gz"
    )


@pytest.fixture(scope="session")
def pgs999999() -> pathlib.Path:
    """A fake scoring file made from PGS001229"""
    return (
        pathlib.Path(__file__).parent
        / "data"
        / "structured"
        / "PGS999999_hmPOS_GRCh38.txt.gz"
    )


@pytest.fixture
def n_samples(vcf_path) -> int:
    """Number of samples in tiny 1000G"""
    with pysam.VariantFile(str(vcf_path)) as vcf:
        return len(vcf.header.samples)


@pytest.fixture
def sample_ids(vcf_path) -> list[str]:
    """Number of samples in tiny 1000G"""
    with pysam.VariantFile(str(vcf_path)) as vcf:
        return list(vcf.header.samples)


@pytest.fixture
def positions_PGS003162():
    """500,000 randomly sampled positions from PGS003162"""
    path = pathlib.Path(__file__).parent / "data" / "PGS003162_positions.txt.gz"
    with gzip.open(path, "rt") as f:
        return list(csv.reader(f, delimiter="\t"))


@pytest.fixture
def x_y_variants(pgs001229):
    positions = Scorefiles(pgs001229).get_unique_positions()
    return list(filter(lambda x: x[0] in {"X", "Y"}, positions))


@pytest.fixture
def pgs001229_autosome(pgs001229):
    positions = Scorefiles(pgs001229).get_unique_positions()
    return list(filter(lambda x: x[0] in {str(x) for x in range(1, 23)}, positions))


@pytest.fixture(scope="session")
def target_with_missingness(vcf_path, tmp_path_factory, missing_prob):
    """A VCF where ~30% of samples have randomly missing genotypes in each variant"""
    cache = tmp_path_factory.mktemp("missing")
    missing_vcf = cache / "missing.vcf.gz"
    missing_prob = missing_prob
    random.seed(42)

    new_recs = []
    with pysam.VariantFile(str(vcf_path), mode="r") as vcf:
        header = vcf.header.copy()
        for rec in vcf.fetch():
            gts = [s["GT"] for s in rec.samples.values()]
            missing_gts = [
                None if random.random() < missing_prob else tup for tup in gts
            ]

            for i, sample in enumerate(rec.samples.values()):
                sample["GT"] = missing_gts[i]

            new_recs.append(rec)

    with pysam.VariantFile(str(missing_vcf), mode="w", header=header) as vcf_out:
        for rec in new_recs:
            vcf_out.write(rec)

    pysam.bcftools.index(str(missing_vcf))

    return missing_vcf


@pytest.fixture(scope="session")
def missing_prob():
    return 0.3


@pytest.fixture
def missing_bgen():
    """This file was recoded from the missing_vcf VCF using plink2 (chrom 1 - 22)"""
    return pathlib.Path(__file__).parent / "data" / "bgen" / "missing" / "missing.bgen"


@pytest.fixture
def missing_bgen_sample():
    """This file was recoded from the missing_vcf VCF using plink2 (chrom 1 - 22)"""
    return (
        pathlib.Path(__file__).parent / "data" / "bgen" / "missing" / "missing.sample"
    )


@pytest.fixture
def gt_array_with_missingness(
    tmp_path_factory, target_with_missingness, test_positions, sampleset_name
):
    cache_dir = tmp_path_factory.mktemp("cache")
    TargetGenome(
        target_path=target_with_missingness,
        cache_dir=cache_dir,
        sampleset=sampleset_name,
    ).cache_variants(positions=test_positions)

    # test parquet files are written
    store = zarr.storage.LocalStore(cache_dir / "gts")
    group = zarr.open_group(
        store=store, path=f"{sampleset_name}/{target_with_missingness.name}"
    )
    array = group["genotypes"]

    return np.array(array[: len(test_positions), :, :])


@pytest.fixture
def effect_allele_idx(gt_array_with_missingness):
    """Randomly generated effect allele indices with the same shape as
    gt_array_with_missingness."""
    rng = np.random.default_rng(42)
    return rng.integers(low=0, high=2, size=gt_array_with_missingness.shape[0])


@pytest.fixture
def genotype_missing_masked(gt_array_with_missingness):
    """An array where sentinel values are replaced with np.nan"""
    masked_array = gt_array_with_missingness.copy()
    masked_array = masked_array.astype(np.float64)
    masked_array[masked_array == MISSING_GENOTYPE_SENTINEL_VALUE] = np.nan
    return masked_array


@pytest.fixture(scope="function")
def missing_dosage_array(genotype_missing_masked, effect_allele_idx, tmp_path):
    missing_array = da.from_array(genotype_missing_masked)
    dosage = calculate_effect_allele_dosage(
        genotype_array=missing_array, effect_idx=effect_allele_idx
    )
    dosage_zarr_path = tmp_path / "dosage"
    dosage.to_zarr(dosage_zarr_path)

    return zarr.open_array(dosage_zarr_path)
