from pathlib import Path

import pytest

from pgscatalog.validate.cli.validate_cli import run
from unittest.mock import patch


def assert_valid_output(captured_stdout):
    assert_invalid_output(captured_stdout, 0, 0)


def assert_invalid_output(captured_stdout, expected_n_errors: int, expected_n_warnings: int = 0):
    stdout_list = [line for line in str(captured_stdout.out).splitlines() if not line.startswith('#')]
    errors, warnings = [], []
    for line in stdout_list:
        (warnings if line.startswith('WARNING:') else errors).append(line)
    assert len(errors) == expected_n_errors
    assert len(warnings) == expected_n_warnings


@pytest.mark.parametrize("scoring_file", [
    "data/valid_raw.txt",
    "data/valid_raw.txt.gz",
    "data/valid_formatted.txt.gz"
])
def test_valid(scoring_file, capsys):

    args: list[str] = ["pgscatalog-validate", "-f", str(Path(__file__).parent / scoring_file)]

    with patch("sys.argv", args):
        run()
        captured_stdout = capsys.readouterr()
        assert_valid_output(captured_stdout)


@pytest.mark.parametrize("scoring_file,expected_n_errors", [
    ("data/invalid_raw_5_errors.txt", 5)
])
def test_invalid(scoring_file, expected_n_errors, capsys):

    args: list[str] = ["pgscatalog-validate", "-f", str(Path(__file__).parent / scoring_file)]

    with patch("sys.argv", args):
        run()
        captured_stdout = capsys.readouterr()
        assert_invalid_output(captured_stdout, expected_n_errors)
