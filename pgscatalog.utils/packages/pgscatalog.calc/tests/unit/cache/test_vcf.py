import pathlib

import numpy as np
import pysam
import pytest
import zarr
import zarr.storage

from pgscatalog.calc.lib.cache._vcf import (
    region_sort_key,
    vcf_buffer_variants,
    vcf_get_genotypes,
    vcf_get_sample_list,
)
from pgscatalog.calc.lib.cache.genotypecache import meta_group_path
from pgscatalog.calc.lib.cache.targetgenome import TargetGenome
from pgscatalog.calc.lib.cache.targetvariants import (
    TargetVariants,
)
from pgscatalog.calc.lib.constants import (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    MISSING_GENOTYPE_TUPLE,
    NO_ALLELE_SENTINEL_VALUE,
)

HEMIZYGOUS_NO_CALL = (MISSING_GENOTYPE_SENTINEL_VALUE, NO_ALLELE_SENTINEL_VALUE)

VCF_HEADER = """##fileformat=VCFv4.2
##contig=<ID=1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}
"""


def one_variant(tmp_path, calls: list[str]) -> pysam.VariantRecord:
    """A single-variant VCF written as text, so the GT spellings under test are
    exactly what a writer would have emitted"""
    path = tmp_path / "one.vcf"
    samples = "\t".join(f"sample{i}" for i in range(len(calls)))
    with path.open("w") as f:
        f.write(VCF_HEADER.format(samples=samples))
        f.write("1\t100\t.\tA\tG\t.\t.\t.\tGT\t" + "\t".join(calls) + "\n")

    with pysam.VariantFile(str(path)) as vcf:
        return next(iter(vcf))


def multiallelic_target(
    tmp_path, calls: list[str]
) -> tuple[pathlib.Path, pathlib.Path]:
    """An indexed one-variant target holding `A -> C,G` at 1:100

    Written as text and indexed, so bcftools sees exactly the record a writer would
    have emitted and the filters under test get a chance to drop it.
    """
    path = tmp_path / "multiallelic.vcf"
    samples = "\t".join(f"sample{i}" for i in range(len(calls)))
    with path.open("w") as f:
        f.write(VCF_HEADER.format(samples=samples))
        f.write("1\t100\t.\tA\tC,G\t.\t.\t.\tGT\t" + "\t".join(calls) + "\n")

    compressed = pathlib.Path(pysam.tabix_index(str(path), preset="vcf", force=True))
    return compressed, pathlib.Path(str(compressed) + ".tbi")


@pytest.fixture
def test_positions() -> list[tuple[str, int]]:
    """Some positions in tiny1000G"""
    return [("8", 116938529), ("1", 44465567), ("10", 97611390)]


def test_a_haploid_no_call_is_stored_hemizygous(tmp_path):
    """`.` is a call of ploidy one that says nothing about which allele, so it is
    a hemizygous no-call - the [255, 254] row of ADR 0007's encoding table. It used
    to be flattened to a diploid missing call and imputed at twice the dosage.
    """
    record = one_variant(tmp_path, ["0", "1", "."])

    assert vcf_get_genotypes(record) == [
        (0, NO_ALLELE_SENTINEL_VALUE),
        (1, NO_ALLELE_SENTINEL_VALUE),
        HEMIZYGOUS_NO_CALL,
    ]


def test_a_diploid_no_call_is_stored_diploid(tmp_path):
    """`./.` is two alleles, neither of them known"""
    record = one_variant(tmp_path, ["0/1", "./."])

    assert vcf_get_genotypes(record) == [(0, 1), MISSING_GENOTYPE_TUPLE]


def test_a_variant_nobody_called_has_nothing_to_compare_against(tmp_path):
    """Every call is `.`, so there is no diploid call to contradict it and the
    length is taken at face value"""
    record = one_variant(tmp_path, [".", ".", "."])

    assert vcf_get_genotypes(record) == [HEMIZYGOUS_NO_CALL] * 3


def test_a_dot_where_every_other_call_is_diploid_is_rejected(tmp_path):
    """Reading ploidy from the GT length means trusting the writer, so the trust is
    checked: nothing at this variant is hemizygous, so a `.` has no reading here
    and can only be a writer that dropped an allele. Imputing those samples at half
    the dosage they can carry is the same class of error as the one this fixes, so
    the variant is named rather than guessed at.
    """
    record = one_variant(tmp_path, ["0/1", "1/1", "."])

    with pytest.raises(ValueError) as excinfo:
        vcf_get_genotypes(record)

    message = str(excinfo.value)
    assert "1:100" in message
    assert "Non-conformant VCF" in message


def test_a_dot_is_trusted_where_some_call_is_hemizygous(tmp_path):
    """A mixed-sex chrX variant: a hemizygous call sits right there, so hemizygosity
    is real at this position and a `.` is a plausible reading of it. Refusing this
    would refuse ordinary conformant output for most chrX data.

    The cost is stated in ADR 0007 rather than fixed: a diploid sample's no-call
    miswritten `.` passes here and is imputed at half its ploidy. Separating the two
    needs a declared sample sex, which ADR 0007 refuses.
    """
    record = one_variant(tmp_path, ["0/1", "1", "."])

    assert vcf_get_genotypes(record) == [
        (0, 1),
        (1, NO_ALLELE_SENTINEL_VALUE),
        HEMIZYGOUS_NO_CALL,
    ]


def test_a_partially_missing_call_keeps_its_known_allele(tmp_path):
    """`./1` is diploid and half known: the known allele survives"""
    record = one_variant(tmp_path, ["0/1", "./1"])

    assert vcf_get_genotypes(record) == [
        (0, 1),
        (MISSING_GENOTYPE_SENTINEL_VALUE, 1),
    ]


def test_a_call_wider_than_the_ploidy_axis_is_rejected(tmp_path):
    """Keeping the first two alleles would turn a call the cache can't hold into a
    plausible score (ADR 0007)"""
    record = one_variant(tmp_path, ["0/1/1"])

    with pytest.raises(ValueError, match="Ploidy 3"):
        vcf_get_genotypes(record)


def test_region_sort_key_accepts_non_numeric_chromosomes():
    """int("X") raised, so nothing downstream could run on a non-autosome"""
    regions = [("X", 2, 2), ("10", 1, 1), ("MT", 1, 1), ("2", 9, 9), ("2", 1, 1)]

    assert sorted(regions, key=region_sort_key) == [
        ("2", 1, 1),
        ("2", 9, 9),
        ("10", 1, 1),
        ("MT", 1, 1),
        ("X", 2, 2),
    ]


# Test 1: vcf_get_sample_list
def test_vcf_get_sample_list_reads_header_samples_from_fixture(vcf_path):
    samples = vcf_get_sample_list(vcf_path)
    # Just sanity-check that we got some samples and the type is correct
    assert isinstance(samples, list)
    assert len(samples) == 3202


# Test 2: vcf_buffer_variants <- run_bcftools_view,
#         parse_target_variant <- vcf_get_genotypes
# Test 2.1: normal case
def test_vcf_buffer_variants(
    test_positions, vcf_path, vcf_index_path, tmp_path_factory, sampleset_name
):
    cache_dir = tmp_path_factory.mktemp("cache")
    variants = vcf_buffer_variants(
        position_batch=test_positions,
        target_path=vcf_path,
        index_path=vcf_index_path,
        cache_dir=cache_dir,
        sampleset=sampleset_name,
    )

    # test output type
    assert isinstance(variants, TargetVariants)
    # test genotype properties
    assert variants.genotypes.shape == (3, 3202, 2)
    # test variant_metadata properties
    assert variants.variant_metadata.chr_name == ["1", "8", "10"]
    assert variants.variant_metadata.chr_pos == [44465567, 116938529, 97611390]
    assert variants.variant_metadata.ref == ["G", "G", "A"]
    assert variants.variant_metadata.alts == [["A"], ["C"], ["G"]]
    # test variant_ids
    assert variants.variant_ids == [
        "1:44465567:G:A",
        "8:116938529:G:C",
        "10:97611390:A:G",
    ]


# Test 2.2: missing positions
def test_vcf_buffer_variants_with_missing(
    test_positions,
    vcf_path,
    vcf_index_path,
    tmp_path_factory,
    sampleset_name,
    n_samples,
):
    missing_chrom, missing_pos = "1", 9067156  # definitely missing
    positions = test_positions + [(missing_chrom, missing_pos)]
    cache_dir = tmp_path_factory.mktemp("cache")
    variants = vcf_buffer_variants(
        position_batch=positions,
        target_path=vcf_path,
        index_path=vcf_index_path,
        cache_dir=cache_dir,
        sampleset=sampleset_name,
    )

    assert len(variants) == len(positions)
    missing_index = positions.index((missing_chrom, missing_pos))
    assert variants.variant_metadata.chr_name[missing_index] == "1"
    assert variants.variant_metadata.ref[missing_index] is None
    assert variants.variant_metadata.alts[missing_index] is None

    # test that missing genotypes are filled with MISSING_GENOTYPE_TUPLE
    missing_gts = variants.genotypes[missing_index]
    print(missing_gts.shape)
    assert missing_gts.shape == (n_samples, 2)
    assert missing_gts.dtype == np.uint8
    assert np.all(missing_gts == MISSING_GENOTYPE_TUPLE)


# Test 2.3: with haploid chrX genotypes
def test_vcf_buffer_x_variants(
    x_variants, vcf_path, vcf_index_path, tmp_path_factory, sampleset_name
):
    """Ploidy is read from the length of the GT field (ADR 0007).

    1000G writes a man's non-pseudoautosomal X call as `1`, not `1/1`, so the
    fixture holds 1599 hemizygous calls and 1603 diploid ones at every variant.
    Before this the array was built with `np.fromiter(..., dtype=(np.uint8, 2))`,
    which broadcasts `(1,)` to `[1, 1]` and quietly doubles every man's dosage.
    """
    cache_dir = tmp_path_factory.mktemp("cache")
    variants = vcf_buffer_variants(
        position_batch=x_variants,
        target_path=vcf_path,
        index_path=vcf_index_path,
        cache_dir=cache_dir,
        sampleset=sampleset_name,
    )

    assert set(variants.variant_metadata.chr_name) == {"X"}

    ploidy = np.sum(variants.genotypes != NO_ALLELE_SENTINEL_VALUE, axis=2)
    called = variants.genotypes[ploidy == 1]

    # a hemizygous call keeps its one allele: nothing was duplicated into the
    # slot the sample doesn't have
    assert np.all(np.isin(called[:, 0], [0, 1]))
    assert np.all(called[:, 1] == NO_ALLELE_SENTINEL_VALUE)

    # every variant the target actually holds is called in the same mixture of
    # men and women, which is why ploidy is a property of a call
    genotyped = [ref is not None for ref in variants.variant_metadata.ref]
    assert {tuple(np.bincount(row, minlength=3)[1:]) for row in ploidy[genotyped]} == {
        (1599, 1603)
    }


# Test 2.4: multiallelic records (ADR 0010)
def test_a_multiallelic_record_reaches_the_cache(
    tmp_path, tmp_path_factory, sampleset_name
):
    """`--max-alleles 2` dropped the record before the cache existed, which made
    `--keep_multiallelic` unreachable code (ADR 0010).

    Alleles are positional - a call names the allele it carries by index - so the
    whole ALT column is kept in the order the record declared it, and a `2` still
    means the second alternative allele.
    """
    target, index = multiallelic_target(tmp_path, ["0/0", "0/1", "1/2", "2/2"])

    variants = vcf_buffer_variants(
        position_batch=[("1", 100)],
        target_path=target,
        index_path=index,
        cache_dir=tmp_path_factory.mktemp("cache"),
        sampleset=sampleset_name,
    )

    assert len(variants) == 1
    assert variants.variant_metadata.ref == ["A"]
    assert variants.variant_metadata.alts == [["C", "G"]]
    assert variants.variant_ids == ["1:100:A:C/G"]
    np.testing.assert_array_equal(
        variants.genotypes[0], [[0, 0], [0, 1], [1, 2], [2, 2]]
    )


def test_a_multiallelic_position_is_not_padded_absent(
    tmp_path, tmp_path_factory, sampleset_name
):
    """A filtered record never entered `seen_positions`, so it was padded with
    `ref=None` - and `update_match_table`'s `WHERE target.ref IS NOT NULL` then
    dropped it, telling the user a variant wasn't found in data that held it.

    A position genuinely absent from the target is still padded, which is the row
    this has to be told apart from.
    """
    target, index = multiallelic_target(tmp_path, ["0/1", "1/2"])

    variants = vcf_buffer_variants(
        position_batch=[("1", 100), ("1", 200)],
        target_path=target,
        index_path=index,
        cache_dir=tmp_path_factory.mktemp("cache"),
        sampleset=sampleset_name,
    )

    metadata = variants.variant_metadata
    assessed = metadata.chr_pos.index(100)
    absent = metadata.chr_pos.index(200)

    assert metadata.ref[assessed] == "A"
    assert metadata.alts[assessed] == ["C", "G"]
    assert metadata.ref[absent] is None
    assert metadata.alts[absent] is None


def test_a_multiallelic_call_round_trips_through_the_cache(
    tmp_path, tmp_path_factory, sampleset_name
):
    """A `1/2` call has to survive the write to zarr and the read back out of it:
    dosage counts calls equal to an allele index, so an index that shifts on the way
    through would count a different allele"""
    target, index = multiallelic_target(tmp_path, ["0/0", "1/2"])
    cache_dir = tmp_path_factory.mktemp("cache")

    genome = TargetGenome(
        target_path=target,
        target_index_path=index,
        cache_dir=cache_dir,
        sampleset=sampleset_name,
    )
    genome.cache_variants([("1", 100)])

    store = zarr.storage.LocalStore(cache_dir / "genotypes.zarr")
    group = zarr.open_group(store=store, path=f"{sampleset_name}/ALL/{target.name}")
    meta = zarr.open_group(store=store, path=meta_group_path(group.path))

    alts = [x.decode() if isinstance(x, bytes) else x for x in meta["alts"][...]]
    assert alts == ["C/G"]
    np.testing.assert_array_equal(group["genotypes"][...], [[[0, 0], [1, 2]]])
