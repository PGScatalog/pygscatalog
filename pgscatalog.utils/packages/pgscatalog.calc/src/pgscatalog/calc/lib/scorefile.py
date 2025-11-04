from __future__ import annotations

import logging
import os
import pathlib
from typing import TYPE_CHECKING

import duckdb

if TYPE_CHECKING:
    from .types import Pathish, PathishList

logger = logging.getLogger(__name__)


class Scorefiles:
    """
    One or more scoring files processed with the `pgscatalog-format` program.
    """

    column_types = {
        "chr_name": "VARCHAR",
        "chr_position": "UINTEGER",
        "effect_allele": "VARCHAR",
        "other_allele": "VARCHAR",
        "effect_weight": "DOUBLE",
        "effect_type": "VARCHAR",
        "is_duplicated": "BOOLEAN",
        "accession": "VARCHAR",
        "row_nr": "UINTEGER",
    }

    def __init__(self, paths: Pathish | PathishList):
        score_paths: list[pathlib.Path]
        if isinstance(paths, str | os.PathLike):
            score_paths = [pathlib.Path(paths).resolve(strict=True)]
        else:
            score_paths = [pathlib.Path(p).resolve(strict=True) for p in paths]
        self._score_paths = score_paths

    @property
    def paths(self) -> list[pathlib.Path]:
        """Return the list of scoring file paths."""
        return self._score_paths

    def get_unique_positions(self, chrom: str | None = None) -> list[tuple[str, int]]:
        path_strings: list[str] = [str(pathlib.Path(x).resolve()) for x in self.paths]

        with duckdb.connect() as conn:
            try:
                query = conn.sql(f"""
                    SELECT DISTINCT chr_name, chr_position
                    FROM read_csv({path_strings}, columns={self.column_types})
                    WHERE chr_name IS NOT NULL AND chr_position IS NOT NULL
                """)
            except duckdb.InvalidInputException as e:
                logger.critical(
                    "Couldn't parse the scoring file. "
                    "Did you process the file with pgscatalog-format? "
                    f"Error: {e}"
                )
                raise

            if chrom is not None:
                query = query.filter(f"chr_name = '{chrom}'")

            logger.info(f"Querying unique positions from scoring files: {path_strings}")

            results: list[tuple[str, int]] = query.order(
                "chr_name, chr_position"
            ).fetchall()
            logger.info(f"Found {len(results)} unique positions")
            return results


def load_scoring_files(
    db_path: Pathish, scorefile_paths: PathishList, max_memory_gb: str, threads: int
) -> None:
    """Load scoring files into the score_variant_table

    Parameters
    ----------
    db_path : Pathish
        Path to the DuckDB database file.
    scorefile_paths : PathishList
        A list of Pathish objects to the scoring CSV file(s). Scoring files must be
        in a  structured format as created by pgscatalog-format.
    max_memory_gb : str
        Maximum memory DuckDB is allowed to use (e.g., "4GB").
    threads : int
        Number of threads for DuckDB to use.

    Notes
    -----
    The score_variant_table is created or replaced each time this function is called.

    Effect weights are stored as double precision floating-point numbers
    (np.float64 equivalent). All previous processing (e.g. by pgscatalog.core) treats
    effect weights as strings to prevent precision problems.

    pgscatalog-format can process scoring files from the PGS Catalog or custom scoring
    files.
    """
    path_strings: list[str] = [str(pathlib.Path(x).resolve()) for x in scorefile_paths]

    with duckdb.connect(
        str(db_path),
        config={"max_memory": max_memory_gb, "threads": str(threads)},
    ) as conn:
        logger.info(f"Loading {path_strings} into {db_path=} score_variant_table")
        conn.sql("""
        CREATE OR REPLACE TABLE score_variant_table (
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
        # column types are manually set here, because chr_name often fools auto-detect
        # (1-22... X)
        conn.execute(
            """
        INSERT INTO score_variant_table
        SELECT * FROM read_csv(?, columns={
            'chr_name': 'VARCHAR',
            'chr_position': 'UINTEGER',
            'effect_allele': 'VARCHAR',
            'other_allele': 'VARCHAR',
            'effect_weight': 'DOUBLE',
            'effect_type': 'VARCHAR',
            'is_duplicated': 'BOOLEAN',
            'accession': 'VARCHAR',
            'row_nr': 'UINTEGER'
        });
        """,
            [path_strings],
        )
        logger.info("Finished loading score_variant_table")
