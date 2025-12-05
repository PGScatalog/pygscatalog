import duckdb
import pytest
from decimal import *

from pgscatalog.calc.lib.score._matchlog import (
    add_complement_macro,
    add_match_log,
    get_ok_accessions,
)

from pgscatalog.calc.lib.scorefile import load_scoring_files

from pgscatalog.calc.lib.score._matchvariants import (
    update_match_table,
    _setup_enums,
)

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

@pytest.fixture
def db_path(tmp_path):
    """Path to a temporary DuckDB file used by tests."""
    return tmp_path / "test_matchlog.duckdb"

# test the add_match_log create the table score_log_table
# Test the add_match_log insert the correct data into score_log_table
# Test the add_match_log create the table summary_log_table
def _make_test_conn(tmp_path,pgs000001) -> duckdb.DuckDBPyConnection:
    db_path = tmp_path / "test.duckdb"

    # score_variant_table loading via load_scoring_files functions from Scorefiles
    load_scoring_files(
        db_path=db_path,
        scorefile_paths=[pgs000001],
        threads=1,
        max_memory_gb="4GB",
    )

    conn = duckdb.connect(db_path)

    # minimal targetvariants table matching your SQL in update_match_table
    conn.execute("""
        CREATE TABLE targetvariants (
            geno_index UINTEGER,
            chr_name TEXT,
            chr_pos UINTEGER,
            ref TEXT,
            alts TEXT[],
            filename TEXT,
            sampleset TEXT
        );
    """)

    # One matching target row in the same sampleset
    conn.execute(
        """
        INSERT INTO targetvariants
        (geno_index, chr_name, chr_pos, ref, alts, filename, sampleset)
        VALUES (?, ?, ?, ?, ?, ?, ?);
        """,
        [10, "11", 69516650 , "T", ["C"], "file.bgen", "test"],
    )

    return conn

def test_add_match_log_with_one_mapped_variants(db_path, sampleset_name, pgs000001, tmp_path):
    """add_match_log should create score_log_table and summary_log_table."""
    conn=_make_test_conn(tmp_path,pgs000001)

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
    assert sum_fraction == Decimal('0.01')
    assert sum_match_rate == Decimal('0.01')
    assert sum_is_ok is False

    conn.close()

# Test with more variants
def _make_conn_from_json(
    db_path,
    allele_match_table,
    score_variant_table,
) -> duckdb.DuckDBPyConnection:
    """
    Create a DuckDB DB with:
    - score_variant_table (from JSON)
    - allele_match_table (from JSON)
    - create an empty targetvariants table
    """
    conn = duckdb.connect(str(db_path))
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
    conn.execute(
        """
        CREATE TABLE allele_match_table AS
        SELECT * 
        FROM read_json(?, format='nd');
        """,
        [str(allele_match_table)],
    )

    # targetvariants table is needed, but left empty at here, so only the pre-defined data is used.

    conn.execute(
       """
       CREATE TABLE targetvariants (
           geno_index UINTEGER,
           chr_name TEXT,
           chr_pos UINTEGER,
           ref TEXT,
           alts TEXT[],
           filename TEXT,
           sampleset TEXT
       );
       """
    )

    return conn

def test_add_match_log_with_ndjson(
    db_path,
    allele_match_table,
    score_variant_table,
):
    sampleset = "1000G"

    conn = _make_conn_from_json(db_path, allele_match_table, score_variant_table)

    #add_complement_macro(conn)

    # Test data includes two sample sets, one test data with match_rate=1.0; while 1000G with match_rate=0.93, so we set threshold 0.9 at here.
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

    print(summary_test,summary_1000G)

    assert test_accession == "test"
    assert test_1000G_accession == "PGS001229_hmPOS_GRCh38"
    assert test_match_rate == Decimal('1.0')
    assert test_1000G_match_rate == Decimal('0.93')
    assert test_sum_is_ok is True
    assert test_1000G_sum_is_ok is False
    conn.close()

    # "test" sampleset pass the threshold 0.95
    accessions = get_ok_accessions(db_path, "test")
    assert accessions[0] == "test"

    # Test the get_ok_accessions function returns ValueError for 1000G sampleset
    with pytest.raises(ValueError, match="No scores passed matching threshold"):
        get_ok_accessions(db_path, "1000G")
    
    
