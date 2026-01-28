from pathlib import Path

import pytest

from pgscatalog.validate.lib.validation import ScoringFileValidation


def resolve_file(input_file: str) -> Path:
    return Path(__file__).parent / input_file


@pytest.mark.parametrize("scoring_file,header", [
    ("data/valid_raw.txt", False),
    ("data/valid_raw.xlsx", False),
    ("data/valid_formatted.txt.gz", True),
    ("data/valid_HLA.txt", False)
])
def test_valid(scoring_file, header):
    validation = ScoringFileValidation(resolve_file(scoring_file), header=header)
    assert len(validation.errors) == 0


@pytest.mark.parametrize("scoring_file,header,expected_n_errors", [
    ("data/invalid_raw_5_errors.txt", False, 5),
    ("data/invalid_HLA_3_errors.txt", False, 3)
])
def test_invalid(scoring_file, header, expected_n_errors):

    validation = ScoringFileValidation(resolve_file(scoring_file), header=header)
    assert len(validation.errors) == expected_n_errors


@pytest.mark.parametrize("scoring_file", [
    "data/trailing_space.txt"
])
def test_toggle_strict(scoring_file):
    input_file = resolve_file(scoring_file)

    validation = ScoringFileValidation(input_file, header=False, strict=False)
    assert len(validation.errors) == 0
    assert len(validation.warnings) == 1

    strict_validation = ScoringFileValidation(input_file, header=False, strict=True)
    assert len(strict_validation.errors) == 1
    assert len(strict_validation.warnings) == 0
