from decimal import Decimal, getcontext

import duckdb
import pytest

from pgscatalog.calc.lib.score._matchlog import (
    add_complement_macro,
    add_match_log,
    get_ok_accessions,
)
from pgscatalog.calc.lib.score._matchvariants import (
    _setup_enums,
    update_match_table,
)
from pgscatalog.calc.lib.scorefile import load_scoring_files


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


# Test if the add_match_log create the table score_log_table
# Test if the add_match_log insert the correct data into score_log_table
# Test if the add_match_log create the table summary_log_table
def _make_test_conn(pgs000001) -> duckdb.DuckDBPyConnection:
    # nothing here needs the database on disk
    conn = duckdb.connect(":memory:")

    # score_variant_table loading via load_scoring_files functions from Scorefiles
    load_scoring_files(conn=conn, scorefile_paths=[pgs000001])

    # minimal targetvariants table matching the SQL in update_match_table
    conn.execute("""
        CREATE TABLE targetvariants (
            geno_index UINTEGER,
            chr_name TEXT,
            chr_pos UINTEGER,
            ref TEXT,
            alts TEXT[],
            zarr_group TEXT,
            sampleset TEXT
        );
    """)

    # One matching target row in the same sampleset
    conn.execute(
        """
        INSERT INTO targetvariants
        (geno_index, chr_name, chr_pos, ref, alts, zarr_group, sampleset)
        VALUES (?, ?, ?, ?, ?, ?, ?);
        """,
        [10, "11", 69516650, "T", ["C"], "test/ALL/file.bgen", "test"],
    )

    return conn


def test_add_match_log_with_one_mapped_variants(sampleset_name, pgs000001):
    """add_match_log should create score_log_table and summary_log_table."""
    conn = _make_test_conn(pgs000001)

    # update_match_table from matchvariants to create allele_match_table
    update_match_table(
        conn=conn,
        match_ambiguous=False,
        match_multiallelic=False,
        sampleset="test",
    )

    add_match_log(conn, sampleset=sampleset_name, min_overlap=0.75)

    tables = conn.execute(
        """
            SELECT table_name
            FROM information_schema.tables
            WHERE table_schema = 'main';
            """
    ).fetchall()

    table_names = {table[0] for table in tables}
    assert "score_log_table" in table_names
    assert "summary_log_table" in table_names

    # Check score_log_table content
    row = conn.execute(
        """
        SELECT
            sampleset,
            accession,
            row_nr,
            chr_name,
            chr_position,
            effect_allele,
            other_allele,
            target_ref,
            target_alts,
            is_matched,
            match_type,
            match_summary,
            matched_effect_allele
        FROM score_log_table
        WHERE match_summary = 'matched';
        """
    ).fetchone()
    (
        records_sampleset,
        records_accession,
        records_row_nr,
        records_chr_name,
        records_chr_position,
        records_eff,
        records_oth,
        records_tref,
        records_talts,
        records_is_matched,
        records_match_type,
        records_match_summary,
        records_matched_effect_allele,
    ) = row

    assert records_sampleset == sampleset_name
    assert records_accession == "PGS000001"
    assert records_row_nr == 0
    assert records_chr_name == "11"
    assert records_chr_position == 69516650
    assert records_eff == "T"
    assert records_oth == "C"
    assert records_tref == "T"
    assert records_talts == ["C"]
    assert records_is_matched is True
    assert records_match_type == "REFALT"
    assert records_match_summary == "matched"
    # REFALT branch uses effect_allele directly
    assert records_matched_effect_allele == "T"
    # Check summary_log_table: single matched row, 100% match_rate, OK

    summary = conn.execute(
        """
        SELECT
            sampleset,
            accession,
            match_summary,
            is_ambiguous,
            is_multiallelic,
            count,
            fraction,
            match_rate,
            is_match_rate_ok
        FROM summary_log_table
        WHERE match_summary = 'matched';
        """
    ).fetchone()
    (
        sum_sampleset,
        sum_accession,
        sum_match_summary,
        sum_is_ambiguous,
        sum_is_multiallelic,
        sum_count,
        sum_fraction,
        sum_match_rate,
        sum_is_ok,
    ) = summary

    print(summary)
    getcontext().prec = 2
    assert sum_sampleset == sampleset_name
    assert sum_accession == "PGS000001"
    assert sum_match_summary == "matched"
    assert sum_is_ambiguous is False
    assert sum_is_multiallelic is False
    assert sum_count == 1
    assert sum_fraction == Decimal("0.01")
    assert sum_match_rate == Decimal("0.01")
    assert sum_is_ok is False

    conn.close()


# Test with more variants
def _make_conn_from_json(
    allele_match_table,
    score_variant_table,
) -> duckdb.DuckDBPyConnection:
    """
    Create an in-memory DuckDB database with:
    - score_variant_table (from JSON)
    - allele_match_table (from JSON)
    """
    conn = duckdb.connect(":memory:")
    _setup_enums(conn)

    # Load the pre-defined test data score_variant_table
    conn.execute(
        """
        CREATE TABLE score_variant_table AS
        SELECT *
        FROM read_json(?, format='nd');
        """,
        [str(score_variant_table)],
    )

    # Load the pre-defined test data allele_match_table
    # target alleles are carried on the match table now, and this fixture predates
    # them, so they log as NULL (which is what the old targetvariants join gave)
    conn.execute(
        """
        CREATE TABLE allele_match_table AS
        SELECT *, NULL::TEXT AS target_ref, NULL::TEXT[] AS target_alts
        FROM read_json(?, format='nd');
        """,
        [str(allele_match_table)],
    )

    return conn


def test_add_match_log_with_ndjson(
    allele_match_table,
    score_variant_table,
):
    sampleset = "1000G"

    conn = _make_conn_from_json(allele_match_table, score_variant_table)

    # Test data includes two samplesets: the test data has match_rate=1.0 and
    # 1000G has match_rate=0.93, so the threshold goes at 0.9
    add_match_log(conn, sampleset=sampleset, min_overlap=0.95)

    # check the match_rate from the summary_log_table
    summary_test = conn.execute(
        """
        SELECT
            accession,
            match_rate,
            is_match_rate_ok
        FROM summary_log_table
        WHERE sampleset='test' AND match_summary='matched';
        """
    ).fetchone()
    (
        test_accession,
        test_match_rate,
        test_sum_is_ok,
    ) = summary_test

    summary_1000G = conn.execute(
        """
        SELECT
            accession,
            match_rate,
            is_match_rate_ok
        FROM summary_log_table
        WHERE sampleset='1000G' AND match_summary='matched';
        """
    ).fetchone()
    (
        test_1000G_accession,
        test_1000G_match_rate,
        test_1000G_sum_is_ok,
    ) = summary_1000G

    print(summary_test, summary_1000G)

    assert test_accession == "test"
    assert test_1000G_accession == "PGS001229_hmPOS_GRCh38"
    assert test_match_rate == Decimal("1.0")
    assert test_1000G_match_rate == Decimal("0.93")
    assert test_sum_is_ok is True
    assert test_1000G_sum_is_ok is False

    # "test" sampleset pass the threshold 0.95
    accessions = get_ok_accessions(conn, "test")
    assert accessions[0] == "test"

    # Test the get_ok_accessions function returns ValueError for 1000G sampleset
    with pytest.raises(ValueError, match="No scores passed matching threshold"):
        get_ok_accessions(conn, "1000G")

    conn.close()


def _make_conn_from_rows(n_variants, match_rows) -> duckdb.DuckDBPyConnection:
    """A scoring file of n_variants, where match_rows (row_nr, match_type,
    is_ambiguous, is_multiallelic) matched and every other variant is unmatched."""
    conn = duckdb.connect(":memory:")
    _setup_enums(conn)
    conn.execute(
        """
        CREATE TABLE score_variant_table AS
        SELECT 'test' AS accession, i AS row_nr, '1' AS chr_name,
               i AS chr_position, 'A' AS effect_allele, 'C' AS other_allele
        FROM generate_series(0, ?) AS t(i);
        """,
        [n_variants - 1],
    )
    conn.execute("""
        CREATE TABLE allele_match_table (
            sampleset TEXT,
            accession TEXT,
            row_nr UINTEGER,
            match_result STRUCT(is_matched BOOLEAN,
                                match_type match_type_enum,
                                match_summary match_summary_enum,
                                is_ambiguous BOOLEAN,
                                is_multiallelic BOOLEAN,
                                match_priority UINTEGER),
            target_row_nr UINTEGER,
            zarr_group TEXT,
            target_ref TEXT,
            target_alts TEXT[]
        );
    """)
    conn.executemany(
        """
        INSERT INTO allele_match_table VALUES
        ('test', 'test', ?,
         {'is_matched': TRUE, 'match_type': ?, 'match_summary': 'matched',
          'is_ambiguous': ?, 'is_multiallelic': ?, 'match_priority': 6},
         0, 'test/ALL/file.vcf.gz', 'A', ['C']);
        """,
        match_rows,
    )
    return conn


def test_match_rate_sums_unrounded_fractions():
    """Rounding each subgroup's fraction before summing them can move the match rate
    across min_overlap, so the match rate must sum raw fractions and round once."""
    # 735 of 1000 variants matched, split evenly across three subgroups: each
    # subgroup's 0.245 rounds up to 0.25, which sums to 0.75 instead of 0.735
    flags = [(False, False), (True, False), (False, True)]
    rows = [(row_nr, "REFALT", *flags[row_nr % 3]) for row_nr in range(735)]
    conn = _make_conn_from_rows(1000, rows)

    add_match_log(conn, sampleset="test", min_overlap=0.75)

    match_rate, is_ok = conn.execute(
        "SELECT DISTINCT match_rate, is_match_rate_ok FROM summary_log_table"
    ).fetchone()
    assert is_ok is False
    assert match_rate < Decimal("0.75")
    # match_rate is still displayed at two decimal places
    assert match_rate.as_tuple().exponent == -2

    conn.close()


def test_match_flipped_counts_strand_flips():
    """match_type_enum values are uppercase and duckdb's LIKE is case sensitive."""
    rows = [(0, "REFALT_FLIP", False, False), (1, "ALTREF_FLIP", False, False)]
    conn = _make_conn_from_rows(2, rows)

    add_match_log(conn, sampleset="test", min_overlap=0.75)

    assert conn.execute(
        "SELECT match_flipped FROM summary_log_table WHERE match_summary = 'matched'"
    ).fetchone() == (2,)

    conn.close()
