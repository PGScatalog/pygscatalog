"""A declared sex: parsing it, and demanding it (ADR 0014)"""

import argparse
import gzip
import pathlib

import numpy as np
import pytest
import zarr
import zarr.storage

from pgscatalog.calc import Scorefiles
from pgscatalog.calc.cli.load import check_sex_map_required
from pgscatalog.calc.lib.constants import (
    MISSING_GENOTYPE_SENTINEL_VALUE,
    NO_ALLELE_SENTINEL_VALUE,
)
from pgscatalog.calc.lib.sexcheck import (
    _outcome,
    _shard_fractions,
    check_sex_map_covers,
    read_sex_map,
)


def _write_map(path: pathlib.Path, body: str) -> pathlib.Path:
    path.write_text(body)
    return path


def _write_scorefile(path: pathlib.Path, chroms: list[str]) -> pathlib.Path:
    """A minimal scoring file carrying one variant per named chromosome"""
    with gzip.open(path, "wt") as f:
        f.write(
            "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\t"
            "effect_type\tis_duplicated\taccession\trow_nr\n"
        )
        for row_nr, chrom in enumerate(chroms):
            f.write(
                f"{chrom}\t{1000 + row_nr}\tA\tG\t1.0\tadditive\tFalse\t"
                f"PGS999999\t{row_nr}\n"
            )
    return path


class TestReadSexMap:
    def test_a_map_is_two_columns_and_three_codes(self, tmp_path):
        path = _write_map(
            tmp_path / "sex.tsv",
            "sample_id\tsex\nHG00096\t1\nHG00097\t2\nHG00099\t0\n",
        )

        assert read_sex_map(path) == {"HG00096": 1, "HG00097": 2, "HG00099": 0}

    def test_zero_is_a_legal_value_not_a_warning(self, tmp_path):
        """0 is how a user exempts a sample from the check, so it is ordinary input"""
        path = _write_map(tmp_path / "sex.tsv", "sample_id\tsex\nHG00096\t0\n")

        assert read_sex_map(path) == {"HG00096": 0}

    def test_a_code_outside_the_coding_names_the_file_line_and_value(self, tmp_path):
        path = _write_map(
            tmp_path / "sex.tsv", "sample_id\tsex\nHG00096\t1\nHG00097\t3\n"
        )

        with pytest.raises(ValueError) as excinfo:
            read_sex_map(path)

        message = str(excinfo.value)
        assert str(path) in message
        assert "line 3" in message  # the number the user's editor shows
        assert "'3'" in message

    def test_a_non_numeric_code_is_refused_the_same_way(self, tmp_path):
        path = _write_map(tmp_path / "sex.tsv", "sample_id\tsex\nHG00096\tmale\n")

        with pytest.raises(ValueError, match="line 2"):
            read_sex_map(path)

    def test_the_header_is_required_and_named(self, tmp_path):
        """No header sniffing: a headerless file is refused, not guessed at"""
        path = _write_map(tmp_path / "sex.tsv", "HG00096\t1\nHG00097\t2\n")

        with pytest.raises(ValueError, match="header"):
            read_sex_map(path)

    def test_wrong_column_names_are_refused(self, tmp_path):
        path = _write_map(tmp_path / "sex.tsv", "IID\tsex\nHG00096\t1\n")

        with pytest.raises(ValueError, match="header"):
            read_sex_map(path)

    def test_a_sample_declared_twice_is_refused(self, tmp_path):
        """Two declarations are not a declaration: last-write-wins would silently
        pick one"""
        path = _write_map(
            tmp_path / "sex.tsv", "sample_id\tsex\nHG00096\t1\nHG00096\t2\n"
        )

        with pytest.raises(ValueError, match="declared twice"):
            read_sex_map(path)

    def test_an_empty_map_is_refused(self, tmp_path):
        path = _write_map(tmp_path / "sex.tsv", "sample_id\tsex\n")

        with pytest.raises(ValueError, match="no samples"):
            read_sex_map(path)


class TestCoverage:
    def test_a_genotyped_sample_missing_from_the_map_fails_naming_it(self):
        """A mandatory declaration that is silently incomplete is not a
        declaration"""
        with pytest.raises(ValueError) as excinfo:
            check_sex_map_covers(
                sex_map={"HG00096": 1},
                samples=["HG00096", "HG00097"],
                source="chrX.vcf.gz",
            )

        assert "HG00097" in str(excinfo.value)

    def test_a_map_wider_than_the_shard_is_fine(self):
        """Sample batches see a subset of the cohort, so one cohort-wide map is the
        normal way to run this (ADR 0003)"""
        check_sex_map_covers(
            sex_map={"HG00096": 1, "HG00097": 2, "HG00099": 0},
            samples=["HG00096"],
            source="batch1.vcf.gz",
        )

    def test_a_map_on_the_wrong_id_space_says_so(self):
        """Every sample missing is the signature of a map keyed on other IDs"""
        with pytest.raises(ValueError, match="sample ID space"):
            check_sex_map_covers(
                sex_map={"1_HG00096": 1},
                samples=["HG00096"],
                source="chrX.vcf.gz",
            )


class TestMandate:
    def test_x_variants_require_a_map(self, tmp_path):
        scorefile = _write_scorefile(tmp_path / "PGS999999.txt.gz", ["1", "X"])

        with pytest.raises(argparse.ArgumentTypeError, match="--sex_map is required"):
            check_sex_map_required(Scorefiles([scorefile]), sex_map=None)

    def test_y_variants_require_a_map(self, tmp_path):
        scorefile = _write_scorefile(tmp_path / "PGS999999.txt.gz", ["1", "Y"])

        with pytest.raises(argparse.ArgumentTypeError, match="--sex_map is required"):
            check_sex_map_required(Scorefiles([scorefile]), sex_map=None)

    def test_mt_variants_do_not(self, tmp_path):
        """MT is hemizygous in everyone and says nothing about sex"""
        scorefile = _write_scorefile(tmp_path / "PGS999999.txt.gz", ["1", "MT"])

        check_sex_map_required(Scorefiles([scorefile]), sex_map=None)

    def test_autosome_only_scoring_files_are_unaffected(self, tmp_path):
        scorefile = _write_scorefile(tmp_path / "PGS999999.txt.gz", ["1", "8"])

        check_sex_map_required(Scorefiles([scorefile]), sex_map=None)

    def test_a_supplied_map_satisfies_the_mandate(self, tmp_path):
        scorefile = _write_scorefile(tmp_path / "PGS999999.txt.gz", ["X"])

        check_sex_map_required(Scorefiles([scorefile]), sex_map=tmp_path / "sex.tsv")


class TestOutcome:
    """The verdict for one sample, given its fraction (ADR 0014)"""

    def test_a_man_hemizygous_throughout_is_consistent(self):
        assert _outcome(declared=1, fraction=1.0, n_informative=34) == "consistent"

    def test_a_woman_diploid_throughout_is_consistent(self):
        assert _outcome(declared=2, fraction=0.0, n_informative=34) == "consistent"

    def test_the_majority_decides_rather_than_unanimity(self):
        """Ten informative variants with a majority rule tolerates four bad calls"""
        assert _outcome(declared=1, fraction=0.7, n_informative=10) == "consistent"
        assert _outcome(declared=2, fraction=0.3, n_informative=10) == "consistent"

    def test_a_declaration_the_majority_contradicts_is_a_mismatch(self):
        assert _outcome(declared=2, fraction=1.0, n_informative=34) == "mismatch"
        assert _outcome(declared=1, fraction=0.0, n_informative=34) == "mismatch"

    def test_an_exact_tie_is_a_mismatch_for_either_declaration(self):
        """No majority at 0.5. Half the calls say male and half say female, so the
        sample or the file is broken either way - contradictory evidence rather than
        the absent evidence the floor covers"""
        assert _outcome(declared=1, fraction=0.5, n_informative=10) == "mismatch"
        assert _outcome(declared=2, fraction=0.5, n_informative=10) == "mismatch"

    def test_below_the_floor_nothing_concludes_even_when_the_calls_disagree(self):
        """Not concluding is a warning, not a failure: refusing every cohort whose PGS
        lacks X coverage would make the mandate self-defeating"""
        assert _outcome(declared=2, fraction=1.0, n_informative=9) == "not_concluded"

    def test_a_null_fraction_never_concludes(self):
        assert _outcome(declared=1, fraction=None, n_informative=34) == "not_concluded"

    def test_zero_exempts_a_sample_whatever_the_calls_say(self):
        assert _outcome(declared=0, fraction=1.0, n_informative=34) == "exempt"
        assert _outcome(declared=0, fraction=0.0, n_informative=34) == "exempt"


def _shard(
    calls: list[list[list[int]]], x_variants=None
) -> tuple[list[float | None], list[int]]:
    """Fractions for a hand-built shard of (n_variants, n_samples, 2) calls"""
    gts = np.array(calls, dtype=np.uint8)
    store = zarr.storage.MemoryStore()
    root = zarr.open_group(store=store)
    genotypes = root.create_array("genotypes", shape=gts.shape, dtype=np.uint8)
    genotypes[:] = gts
    samples = [f"sample{i}" for i in range(gts.shape[1])]
    return _shard_fractions(genotypes, samples, x_variants)


HET = [0, 1]
HOM_ALT = [1, 1]
HEMI_ALT = [1, NO_ALLELE_SENTINEL_VALUE]
HEMI_REF = [0, NO_ALLELE_SENTINEL_VALUE]
HEMI_NOCALL = [MISSING_GENOTYPE_SENTINEL_VALUE, NO_ALLELE_SENTINEL_VALUE]
DIPLOID_NOCALL = [MISSING_GENOTYPE_SENTINEL_VALUE, MISSING_GENOTYPE_SENTINEL_VALUE]


class TestShardFractions:
    def test_a_mixed_sex_variant_separates_the_two(self):
        fractions, n_informative = _shard([[HEMI_ALT, HET], [HEMI_REF, HOM_ALT]])

        assert n_informative == [2, 2]
        assert fractions == [1.0, 0.0]

    def test_a_variant_diploid_in_everyone_is_not_informative(self):
        """Either pseudoautosomal or a cohort with nobody hemizygous in it, and
        nothing tells those apart without PAR coordinates (ADR 0007 declined them)"""
        fractions, n_informative = _shard([[HET, HET], [HEMI_ALT, HOM_ALT]])

        assert n_informative == [1, 1]
        assert fractions == [1.0, 0.0]

    def test_an_all_female_shard_cannot_be_checked_and_says_so(self):
        """No hemizygous call anywhere means no informative variant anywhere. The
        fraction is NULL, never 0.0: 0/0 read as diploid would pass every woman
        silently and flag every man, on a check that never ran (ADR 0014)"""
        fractions, n_informative = _shard([[HET, HOM_ALT], [HET, HET]])

        assert n_informative == [0, 0]
        assert fractions == [None, None]

    def test_a_hemizygous_no_call_makes_a_variant_informative(self):
        """The file declared ploidy one there, which is the signal being read"""
        _, n_informative = _shard([[HEMI_NOCALL, HET]])

        assert n_informative == [0, 1]

    def test_a_hemizygous_no_call_is_in_neither_half_of_its_own_fraction(self):
        """Informative for the variant, absent from the sample's own ratio: the
        denominator is its *called* genotypes"""
        fractions, _ = _shard([[HEMI_NOCALL, HET], [HEMI_ALT, HET]])

        # one called genotype, hemizygous, so 1/1 rather than 1/2
        assert fractions[0] == 1.0

    def test_a_sample_with_no_called_genotype_gets_a_null_fraction(self):
        fractions, n_informative = _shard([[HEMI_ALT, DIPLOID_NOCALL]])

        assert n_informative == [1, 0]
        assert fractions == [1.0, None]


def test_informative_floor_counts_each_samples_called_variants():
    fractions, counts = _shard([[HEMI_ALT, HET]] + [[HEMI_ALT, DIPLOID_NOCALL]] * 9)
    assert counts == [10, 1]
    assert _outcome(1, fractions[1], counts[1]) == "not_concluded"


def test_whole_genome_check_excludes_other_chromosomes():
    fractions, counts = _shard(
        [[HEMI_ALT, HET]] * 10 + [[HEMI_ALT, HEMI_ALT]] * 20,
        x_variants=np.array([True] * 10 + [False] * 20),
    )
    assert fractions == [1.0, 0.0]
    assert counts == [10, 10]
