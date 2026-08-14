"""A declared sex: parsing it, and demanding it (ADR 0014)"""

import argparse
import gzip
import pathlib

import pytest

from pgscatalog.calc import Scorefiles
from pgscatalog.calc.cli.load import check_sex_map_required
from pgscatalog.calc.lib.sexcheck import check_sex_map_covers, read_sex_map


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
