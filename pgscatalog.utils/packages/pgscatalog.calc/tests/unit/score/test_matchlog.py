from decimal import Decimal

import duckdb
import pytest

from pgscatalog.calc.lib.score._matchlog import (
    add_match_log,
    add_reverse_complement_macro,
    get_ok_accessions,
)
from pgscatalog.calc.lib.score._matchvariants import (
    _setup_enums,
    reverse_complement,
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
        ("ATGC", "GCAT"),
        ("AT", "AT"),
        ("A" * 200, "T" * 200),
        ("atgc", "cgta"),
    ],
    ids=[
        "A reverse complements to T",
        "T reverse complements to A",
        "C reverse complements to G",
        "G reverse complements to C",
        # a plain complement gives TACG here: the reversal is the whole bug
        "multiple mixed nts are reversed as well as complemented",
        "a self-reverse-complementary sequence is unchanged",
        "very long sequence",
        # the macro used to upper() and the Python function never has. Neither does
        # now: a lowercase allele can't reach a flipped match type, so upper-casing
        # could only print a sequence the matcher never compared
        "a lowercase allele is left alone, exactly as the matcher leaves it",
    ],
)
def test_reverse_complement_macro(input, result):
    """The SQL macro and the Python function have to give the same answer.

    matched_effect_allele is derived with the macro, so if they disagree the match
    log names a different allele from the match that produced it (docs/adr/0011).
    """
    with duckdb.connect(":memory:") as conn:
        add_reverse_complement_macro(conn)
        assert conn.sql(f"SELECT reverse_complement('{input}')").fetchone() == (result,)
    assert reverse_complement(input) == result


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
        match_flipped=True,
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
                                exclusion_reason exclusion_reason_enum,
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
          'exclusion_reason': NULL,
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


# the summary log is where "217 variants excluded" has to become counts a user can
# act on: raise --keep_multiallelic, raise --keep_ambiguous, drop
# --ignore_strand_flips, or look at the scoring file (docs/adr/0010)

_SHARD = "test/X/x.vcf.gz"
_REASON_SCORE_ROWS = [
    # two score variants resolving to the reference allele of target row 5
    ("X", 100, "A", "C", 0.1, "additive", False, "PGS000001", 0),
    ("X", 100, "A", "C", 0.2, "additive", False, "PGS000001", 1),
    # a palindromic pairing, then a target carrying two alternative alleles
    ("X", 200, "A", "T", 0.3, "additive", False, "PGS000001", 2),
    ("X", 300, "A", "C", 0.4, "additive", False, "PGS000001", 3),
    # a clean match, and a score variant that matches nothing
    ("X", 400, "G", "T", 0.5, "additive", False, "PGS000001", 4),
    ("X", 500, "A", "G", 0.6, "additive", False, "PGS000001", 5),
    # matches target row 10 only when read complemented: revcomp(T/G) is A/C
    ("X", 600, "T", "G", 0.7, "additive", False, "PGS000001", 6),
]
_REASON_TARGET_ROWS = [
    (5, "X", 100, "A", ["C"], _SHARD, "test"),
    (6, "X", 200, "A", ["T"], _SHARD, "test"),
    (7, "X", 300, "A", ["C", "G"], _SHARD, "test"),
    (8, "X", 400, "G", ["T"], _SHARD, "test"),
    (9, "X", 500, "G", ["T"], _SHARD, "test"),
    (10, "X", 600, "A", ["C"], _SHARD, "test"),
]


def _summary_of_every_reason(
    *, match_multiallelic: bool, match_flipped: bool = True
) -> list[tuple]:
    """Match a scoring file that hits all four exclusion reasons, and summarise it"""
    conn = duckdb.connect(":memory:")
    conn.execute("""
        CREATE TABLE score_variant_table (
            chr_name VARCHAR,
            chr_position UINTEGER,
            effect_allele VARCHAR,
            other_allele VARCHAR,
            effect_weight DOUBLE,
            effect_type VARCHAR,
            is_duplicated BOOLEAN,
            accession VARCHAR,
            row_nr UINTEGER,
            PRIMARY KEY(row_nr, accession)
        );
    """)
    conn.executemany(
        "INSERT INTO score_variant_table VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);",
        _REASON_SCORE_ROWS,
    )
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
    conn.executemany(
        "INSERT INTO targetvariants VALUES (?, ?, ?, ?, ?, ?, ?);",
        _REASON_TARGET_ROWS,
    )
    update_match_table(
        conn=conn,
        match_ambiguous=False,
        match_multiallelic=match_multiallelic,
        match_flipped=match_flipped,
        sampleset="test",
    )
    add_match_log(conn, sampleset="test", min_overlap=0.1)

    # cast to sort alphabetically: both columns are enums, so their natural order is
    # the order they were declared in
    summary = conn.execute("""
        SELECT match_summary, exclusion_reason, is_multiallelic, count
        FROM summary_log_table
        ORDER BY match_summary::VARCHAR, exclusion_reason::VARCHAR, is_multiallelic;
    """).fetchall()
    conn.close()
    return summary


def test_summary_log_splits_exclusions_by_reason():
    """One 'excluded' count is several decisions, and each has a different fix"""
    assert _summary_of_every_reason(match_multiallelic=False) == [
        ("excluded", "ambiguous", False, 1),
        ("excluded", "collision", False, 2),
        ("excluded", "multiallelic", True, 1),
        ("matched", None, False, 2),
        ("unmatched", None, False, 1),
    ]


def test_ignoring_strand_flips_adds_a_reason_of_its_own():
    """A flipped match becomes an excluded row that names the flag that dropped it,
    not an unmatched row: the variant was found, and policy discarded it
    (docs/adr/0011)."""
    assert _summary_of_every_reason(match_multiallelic=False, match_flipped=False) == [
        ("excluded", "ambiguous", False, 1),
        ("excluded", "collision", False, 2),
        ("excluded", "multiallelic", True, 1),
        ("excluded", "strand_flip", False, 1),
        ("matched", None, False, 1),
        ("unmatched", None, False, 1),
    ]


def test_keeping_multiallelics_leaves_the_other_reasons_alone():
    """Raising the flag empties one reason's count without touching the others, and
    the variant stays flagged multiallelic now that it's matched: the flag describes
    the target variant, the reason describes the decision."""
    assert _summary_of_every_reason(match_multiallelic=True) == [
        ("excluded", "ambiguous", False, 1),
        ("excluded", "collision", False, 2),
        ("matched", None, False, 2),
        ("matched", None, True, 1),
        ("unmatched", None, False, 1),
    ]


def test_match_flipped_counts_strand_flips():
    """match_type_enum values are uppercase and duckdb's LIKE is case sensitive."""
    rows = [(0, "REFALT_FLIP", False, False), (1, "ALTREF_FLIP", False, False)]
    conn = _make_conn_from_rows(2, rows)

    add_match_log(conn, sampleset="test", min_overlap=0.75)

    assert conn.execute(
        "SELECT match_flipped FROM summary_log_table WHERE match_summary = 'matched'"
    ).fetchone() == (2,)

    conn.close()


# matched_effect_allele goes through the SQL macro, and matching goes through the
# Python function. These run the whole way from a scoring file row to a log row, so a
# macro that disagreed with the matcher would show up here rather than in neither
# test (docs/adr/0011).

_FLIP_SCORE_ROWS = [
    # a SNP read on the other strand: revcomp(A) = T is the reference allele
    ("1", 100, "A", "C", 0.1, "additive", False, "PGS000001", 0),
    # the same score variant where revcomp(A) = T is the *alternative* allele. Which
    # target allele it resolved to differs, the score allele that was read does not
    ("1", 200, "A", "C", 0.2, "additive", False, "PGS000001", 1),
    # a genuine strand flip of a multi-base target: revcomp(AC) = GT, revcomp(CC) = GG
    ("1", 300, "AC", "CC", 0.3, "additive", False, "PGS000001", 2),
    # the plain complement of the same score variant, complement(AC) = TG. This used
    # to match REFALT_FLIP, against a target that is not what the other strand carries
    ("1", 400, "AC", "CC", 0.4, "additive", False, "PGS000001", 3),
]
_FLIP_TARGET_ROWS = [
    (5, "1", 100, "T", ["G"], _SHARD, "test"),
    (6, "1", 200, "G", ["T"], _SHARD, "test"),
    (7, "1", 300, "GT", ["GG"], _SHARD, "test"),
    (8, "1", 400, "TG", ["GG"], _SHARD, "test"),
]


def test_matched_effect_allele_names_the_allele_the_match_resolved_to():
    conn = duckdb.connect(":memory:")
    conn.execute("""
        CREATE TABLE score_variant_table (
            chr_name VARCHAR,
            chr_position UINTEGER,
            effect_allele VARCHAR,
            other_allele VARCHAR,
            effect_weight DOUBLE,
            effect_type VARCHAR,
            is_duplicated BOOLEAN,
            accession VARCHAR,
            row_nr UINTEGER,
            PRIMARY KEY(row_nr, accession)
        );
    """)
    conn.executemany(
        "INSERT INTO score_variant_table VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);",
        _FLIP_SCORE_ROWS,
    )
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
    conn.executemany(
        "INSERT INTO targetvariants VALUES (?, ?, ?, ?, ?, ?, ?);",
        _FLIP_TARGET_ROWS,
    )
    update_match_table(
        conn=conn,
        match_ambiguous=False,
        match_multiallelic=False,
        match_flipped=True,
        sampleset="test",
    )
    add_match_log(conn, sampleset="test", min_overlap=0.1)

    log = conn.execute("""
        SELECT row_nr, match_type, target_ref, target_alts, matched_effect_allele
        FROM score_log_table ORDER BY row_nr;
    """).fetchall()
    conn.close()

    assert log == [
        # the effect allele resolved to the reference allele of the target
        (0, "REFALT_FLIP", "T", ["G"], "T"),
        # ...and here to its alternative allele. The macro used to read the *other*
        # allele in this branch and report G, the allele this variant didn't match
        (1, "ALTREF_FLIP", "G", ["T"], "T"),
        (2, "REFALT_FLIP", "GT", ["GG"], "GT"),
        # a plain complement is not a strand flip: the target variant is still logged,
        # because the position matched, but no allele of it was resolved to
        (3, "NO_MATCH", "TG", ["GG"], None),
    ]
