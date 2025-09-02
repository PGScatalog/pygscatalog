from pathlib import Path

import pytest

from pgscatalog.validate.cli.validate_cli import run
from unittest.mock import patch


def assert_valid_output(captured_stdout):
    assert len(str(captured_stdout).split("\n")) == 1


@pytest.mark.parametrize("scoring_file", [
    "data/valid_raw.txt",
    "data/valid_raw.txt.gz"
])
def test_valid(scoring_file, capsys):

    args: list[str] = ["pgscatalog-validate", "-f", str(Path(__file__).parent / scoring_file)]

    with patch("sys.argv", args):
        run()
        captured_stdout = capsys.readouterr()
        assert_valid_output(captured_stdout)
