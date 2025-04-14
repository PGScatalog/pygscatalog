import subprocess
import pytest


@pytest.mark.parametrize(
    "cli_name",
    [
        "pgscatalog-download",
        "pgscatalog-combine",
        "pgscatalog-format",
        "pgscatalog-relabel",
        "pgscatalog-match",
        "pgscatalog-matchmerge",
        "pgscatalog-intersect",
        "pgscatalog-aggregate",
        "pgscatalog-ancestry-adjust",
    ],
)
def test_clis(cli_name):
    """Test CLI applications are OK in this package"""
    result = subprocess.run(
        [cli_name, "--help"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    assert result.returncode == 0
    assert cli_name in result.stdout.split(" ")[1]  # usage: pgscatalog-download ...
