import sys

import pytest

from pgscatalog.calc.cli.main import main


@pytest.mark.parametrize(
    "argv, expected_text",
    [
        (["pgsc_calc", "--help"], ["usage: pgsc_calc [-h]"]),
        (["pgsc_calc", "load", "--help"], ["usage: pgsc_calc load [-h]"]),
        (["pgsc_calc", "score", "--help"], ["usage: pgsc_calc score [-h]"]),
    ],
)
def test_cli_help(capsys, monkeypatch, argv, expected_text):
    """Ensure each CLI help command runs successfully and prints expected text."""
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        main()

    assert excinfo.value.code == 0

    out, _ = capsys.readouterr()
    for text in expected_text:
        assert text in out
