from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import duckdb

if TYPE_CHECKING:
    from pgscatalog.calc.lib.types import Pathish

logger = logging.getLogger(__name__)


def add_match_log(
    conn: duckdb.DuckDBPyConnection,
    sampleset: str,
    min_overlap: float = 0.75,
) -> None:
    """Describe matching outcome for each scoring file variant
    These logs can get huge - so create a view
    """
    logger.info("Adding score_log_table to database")
    add_complement_macro(conn)

    conn.execute("""
    CREATE TABLE IF NOT EXISTS score_log_table(
        sampleset VARCHAR,
        accession VARCHAR,
        row_nr UINTEGER,
        chr_name VARCHAR,
        chr_position UINTEGER,
        effect_allele VARCHAR NOT NULL,
        matched_effect_allele VARCHAR
            GENERATED ALWAYS AS (
            CASE match_type
                WHEN 'REFALT'       THEN effect_allele
                WHEN 'REF_NO_OA'    THEN effect_allele
                WHEN 'ALTREF'       THEN effect_allele
                WHEN 'ALT_NO_OA'    THEN effect_allele
                WHEN 'REFALT_FLIP'  THEN complement(effect_allele)
                WHEN 'ALTREF_FLIP'  THEN complement(other_allele)
                ELSE NULL
            END
        ),
        other_allele VARCHAR,
        target_ref VARCHAR,
        target_alts VARCHAR[],
        is_matched BOOLEAN NOT NULL,
        match_type match_type_enum NOT NULL,
        match_summary match_summary_enum NOT NULL,
        is_ambiguous BOOLEAN,
        is_multiallelic BOOLEAN,
        target_row_nr UINTEGER,
        filename VARCHAR,
        -- each variant in a scoring file can only match once
        PRIMARY KEY(sampleset, accession, row_nr)
    );
    """)

    # TODO @nebfield: test implementation with multiple samplesets  # noqa: TD003
    conn.execute(
        """
     INSERT INTO score_log_table
     SELECT
         COALESCE(matches.sampleset, ?) AS sampleset,
         score.accession,
         score.row_nr,
         score.chr_name,
         score.chr_position,
         score.effect_allele,
         score.other_allele,
         target.ref AS target_ref,
         target.alts AS target_alts,
         COALESCE(matches.match_result.is_matched, FALSE) AS is_matched,
         COALESCE(matches.match_result.match_type, 'NO_MATCH') AS match_type,
         COALESCE(matches.match_result.match_summary, 'unmatched') AS match_summary,
         matches.match_result.is_ambiguous,
         matches.match_result.is_multiallelic,
         matches.target_row_nr,
         matches.filename
     FROM score_variant_table score
              LEFT JOIN allele_match_table matches
                        ON score.accession = matches.accession
                            AND score.row_nr = matches.row_nr
              LEFT JOIN targetvariants target
                        ON target.geno_index = matches.target_row_nr
                            AND target.filename = matches.filename
     ORDER BY sampleset, score.accession, score.row_nr;
    """,
        [sampleset],
    )

    # now do some calculations to create a nicer summary table
    logger.info("Creating temporary tables to summarise match logs")
    conn.sql(
        """
        CREATE TEMP TABLE temp_summary_log_table AS
        WITH match_counts AS (
            SELECT sampleset,
                   accession,
                   match_summary,
                   is_ambiguous,
                   is_multiallelic,
                   COUNT(match_summary) AS count,
                   -- count any match types which flipped the effect allele
                   COUNT(*) FILTER (WHERE match_type LIKE '%flip') AS match_flipped
            FROM score_log_table
            GROUP BY sampleset,
                     accession,
                     match_summary,
                     is_ambiguous,
                     is_multiallelic
        )
        SELECT *,
               CAST(count / sum(count) OVER (PARTITION BY sampleset, accession)
                   AS DECIMAL(5, 2)) AS fraction
        FROM match_counts;
        """
    )
    conn.execute(
        """
        CREATE TEMP TABLE temp_match_rates_ok AS
        SELECT
            sampleset,
            accession,
            SUM(fraction) AS match_rate,
            SUM(fraction) > ? AS is_match_rate_ok
        FROM temp_summary_log_table
        WHERE match_summary = 'matched'
        GROUP BY sampleset, accession;
        """,
        [min_overlap],
    )

    logger.info("Creating summary_log_table if it doesn't exist")
    conn.sql(
        """
        CREATE TABLE summary_log_table AS
        SELECT
            sampleset,
            accession,
            match_summary,
            is_ambiguous,
            is_multiallelic,
            match_flipped,
            count,
            fraction,
            match_rate,
            is_match_rate_ok
        FROM temp_summary_log_table
        JOIN temp_match_rates_ok
        USING (sampleset, accession);
        """
    )


def get_ok_accessions(db_path: Pathish, sampleset: str) -> list[str]:
    """Get accessions that have passed the match rate threshold.

    Raises a ValueError if no accessions pass the match rate threshold.
    """
    with duckdb.connect(str(db_path), read_only=True) as conn:
        accessions: list[tuple[str, bool, float]] = conn.execute(
            """
            SELECT DISTINCT accession, is_match_rate_ok, match_rate
            FROM summary_log_table
            WHERE sampleset = ?
        """,
            [sampleset],
        ).fetchall()
    passed: list[tuple[str, bool, float]] = list(
        filter(lambda x: x[1] is True, accessions)
    )
    failed: list[tuple[str, bool, float]] = list(
        filter(lambda x: x[1] is False, accessions)
    )

    for accession, _, match_rate in passed:
        logger.info(f"Match rate OK for {accession=} {match_rate=:.2%}")

    for accession, _, match_rate in failed:
        logger.warning(
            f"Not calculating score for {accession=} because {match_rate=:.2%}"
        )
    if len(passed) == 0:
        raise ValueError("No scores passed matching threshold")

    return [x[0] for x in passed]


def add_complement_macro(conn: duckdb.DuckDBPyConnection) -> None:
    conn.sql("""
    CREATE MACRO IF NOT EXISTS complement(allele) AS (
        translate(upper(allele), 'ATCG', 'TAGC')
    );
    """)
