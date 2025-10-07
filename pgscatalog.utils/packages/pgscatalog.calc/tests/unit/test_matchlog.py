import duckdb
import pytest

from pgscatalog.calc.lib._matchlog import add_complement_macro


@pytest.mark.parametrize(
    "input, result",
    [
        ("A", "T"),
        ("T", "A"),
        ("C", "G"),
        ("G", "C"),
        ("ATGC", "TACG"),
        ("A" * 200, "T" * 200),
    ],
    ids=[
        "A complement is T",
        "T complement is A",
        "C complement is G",
        "G complement is C",
        "multiple mixed nts",
        "very long sequence",
    ],
)
def test_complement_macro(input, result):
    with duckdb.connect(":memory:") as conn:
        add_complement_macro(conn)
        assert conn.sql(f"SELECT complement('{input}')").fetchone() == (result,)
