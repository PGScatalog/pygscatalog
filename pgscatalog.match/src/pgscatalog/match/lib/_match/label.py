""" This module labels match candidates with various flag columns, with the aim of
producing a final set of best match candidates (one maximum for each variant).

These operations are all quite data-framey and involve calculating a horizontal
maximum across boolean columns to determine if a variant should be given a true value
in an exclude column.
"""

import logging

import polars as pl

from .preprocess import complement_valid_alleles

logger = logging.getLogger(__name__)


def label_matches(
    df: pl.LazyFrame,
    keep_first_match,
    remove_ambiguous,
    remove_multiallelic,
    skip_flip,
    filter_IDs,
) -> pl.LazyFrame:
    """Label match candidates with additional metadata. Column definitions:

    - match_candidate: All input variants that were returned from match.get_all_matches() (always True in this function)
    - best_match: True if row is the best possible match type (refalt > altref > ...)
    - duplicate: True if more than one best match exists for the same accession and ID
    - ambiguous: True if ambiguous
    """
    labelled = (
        df.with_columns(
            exclude=pl.lit(False)
        )  # set up dummy exclude column for _label_*
        .pipe(_label_best_match)
        .pipe(_label_duplicate_best_match)
        .pipe(_label_duplicate_id, keep_first_match)
        .pipe(_label_biallelic_ambiguous, remove_ambiguous)
        .pipe(_label_multiallelic, remove_multiallelic)
        .pipe(_label_flips, skip_flip)
        .pipe(_label_filter, filter_IDs)
        .with_columns(match_candidate=pl.lit(True))
    )

    return labelled.pipe(_encode_match_priority)


def _encode_match_priority(df: pl.LazyFrame) -> pl.LazyFrame:
    """Encode a new column called match status containing matched, unmatched, excluded, and not_best"""
    return (
        df.with_columns(
            # set false best match to not_best
            match_priority=pl.col("best_match").map_elements(
                lambda x: {None: 0, True: 1, False: 3}[x]
            )
        )
        .with_columns(
            excluded_match_priority=pl.col("exclude").map_elements(
                lambda x: {None: 0, True: 2, False: 0}[x]
            )
        )
        .with_columns(
            max=pl.max_horizontal("match_priority", "excluded_match_priority")
        )
        .with_columns(
            match_status=pl.col("max")
            .map_elements(
                lambda x: {0: "unmatched", 1: "matched", 2: "excluded", 3: "not_best"}[
                    x
                ]
            )
            .cast(pl.Categorical)
        )
        .drop(["max", "excluded_match_priority", "match_priority"])
    )


def _label_best_match(df: pl.LazyFrame) -> pl.LazyFrame:
    """Best matches have the lowest match priority type. Find the best matches and label them."""
    logger.debug("Labelling best match type (refalt > altref > ...)")
    match_priority = {
        "refalt": 0,
        "altref": 1,
        "refalt_flip": 2,
        "altref_flip": 3,
        "no_oa_ref": 4,
        "no_oa_alt": 5,
        "no_oa_ref_flip": 6,
        "no_oa_alt_flip": 7,
    }
    match_priority_rev = {v: k for k, v in match_priority.items()}

    # use a groupby aggregation to guarantee the number of rows stays the same
    # rows were being lost using an anti join + reduce approach
    prioritised: pl.LazyFrame = (
        df.with_columns(match_priority=pl.col("match_type").replace(match_priority))
        .with_columns(
            best_match_type=pl.col("match_priority")
            .min()
            .over(["accession", "row_nr"])
            .replace(match_priority_rev)
        )
        .with_columns(
            best_match=pl.when(pl.col("best_match_type") == pl.col("match_type"))
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
        )
    )

    return prioritised.drop(["match_priority", "best_match_type"])


def _label_duplicate_best_match(df: pl.LazyFrame) -> pl.LazyFrame:
    """A scoring file row_nr in an accession group can be duplicated if a target position has different REF, e.g.:

    ┌────────┬────────────────────────┬────────────┬────────────────┬─────┬────────────┐
    │ row_nr ┆ accession              ┆ match_type ┆ ID             ┆ REF ┆ best_match │
    │ ---    ┆ ---                    ┆ ---        ┆ ---            ┆ --- ┆ ---        │
    │ i64    ┆ cat                    ┆ str        ┆ cat            ┆ str ┆ bool       │
    ╞════════╪════════════════════════╪════════════╪════════════════╪═════╪════════════╡
    │ 38557  ┆ PGS000012_hmPOS_GRCh37 ┆ no_oa_alt  ┆ 3:29588979:A:G ┆ A   ┆ true       │
    ├╌╌╌╌╌╌╌╌┼╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌┼╌╌╌╌╌╌╌╌╌╌╌╌┤
    │ 38557  ┆ PGS000012_hmPOS_GRCh37 ┆ no_oa_alt  ┆ 3:29588979:T:G ┆ T   ┆ true       │
    └────────┴────────────────────────┴────────────┴────────────────┴─────┴────────────┘

    Label the first row with best_match = true, and duplicate rows with best_match = false
    """
    logger.debug(
        "Labelling duplicated best match: keeping first instance as best_match = True"
    )
    labelled: pl.LazyFrame = (
        df.with_columns(
            count=pl.col("best_match").len().over(["accession", "row_nr", "best_match"])
        )
        .with_columns(
            duplicate_best_match=pl.when((pl.col("count") > 1) & (pl.col("best_match")))
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
        )
        .drop("count")
        .with_row_index(
            name="temp_row_nr"
        )  # add temporary row count to get first variant
        .with_columns(
            best_match=pl.when(
                pl.col("best_match")
                & pl.col("duplicate_best_match")
                & (pl.col("temp_row_nr") > pl.min("temp_row_nr")).over(
                    ["accession", "row_nr"]
                )
            )
            .then(False)  # reset best match flag for duplicates
            .otherwise(pl.col("best_match"))  # just keep value from existing column
        )
    )

    return labelled


def _label_duplicate_id(df: pl.LazyFrame, keep_first_match: bool) -> pl.LazyFrame:
    """Label best match duplicates made when the scoring file is remapped to a different genome build

    ┌─────────┬────────────────────────┬─────────────┬────────────────┬─────┬────────────┐
    │ row_nr  ┆ accession              ┆ match_type  ┆ ID             ┆ REF ┆ best_match │
    │ ---     ┆ ---                    ┆ ---         ┆ ---            ┆ --- ┆ ---        │
    │ i64     ┆ cat                    ┆ str         ┆ cat            ┆ str ┆ bool       │
    ╞═════════╪════════════════════════╪═════════════╪════════════════╪═════╪════════════╡
    │ 1194115 ┆ PGS002244_hmPOS_GRCh37 ┆ altref      ┆ 3:50924580:C:A ┆ C   ┆ true       │
    ├╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌┼╌╌╌╌╌╌╌╌╌╌╌╌┤
    │ 1194132 ┆ PGS002244_hmPOS_GRCh37 ┆ refalt_flip ┆ 3:50924580:C:A ┆ C   ┆ true       │
    └─────────┴────────────────────────┴─────────────┴────────────────┴─────┴────────────┘

    refalt > altref > ... prioritisation doesn't fix this problem because row_nr is different (duplicated by remapping)
    """
    logger.debug(
        "Labelling multiple scoring file lines (accession/row_nr) that best_match to the same variant"
    )

    # the window in .over() starts with accession + ID
    # best_match is added to not count: same row_nr, different match_type (_label_best_match)
    # duplicate_best_match is added to not count: same row_nr, same match_type (_label_duplicate_best_match)
    duplicates = df.with_columns(
        count=pl.col("ID")
        .len()
        .over(["accession", "ID", "best_match", "duplicate_best_match"])
    ).with_columns(
        duplicate_ID=pl.when((pl.col("count") > 1) & (pl.col("best_match")))
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
    )

    if keep_first_match:
        logger.debug("Keeping first duplicate, labelling others with exclude flag ")
        # set first duplicate (with the smallest row_nr) to exclude = false
        labelled = duplicates.with_columns(
            exclude_duplicate=pl.when(
                (pl.col("duplicate_ID"))
                & (
                    pl.col("row_nr")
                    != pl.min("row_nr").over(["accession", "ID", "duplicate_ID"])
                )
            )
            .then(True)
            .otherwise(False)
        )
    else:
        logger.debug("Labelling all duplicates with exclude flag")
        labelled = duplicates.with_columns(exclude_duplicate=pl.col("duplicate_ID"))

    # get the horizontal maximum to combine the exclusion columns for each variant
    return (
        labelled.with_columns(max=pl.max_horizontal("exclude", "exclude_duplicate"))
        .drop(["exclude", "exclude_duplicate"])
        .rename({"max": "exclude"})
    )


def _label_biallelic_ambiguous(df: pl.LazyFrame, remove_ambiguous) -> pl.LazyFrame:
    logger.debug("Labelling ambiguous variants")
    ambig = (
        df.with_columns(
            [
                pl.col(
                    [
                        "effect_allele",
                        "other_allele",
                        "REF",
                        "ALT",
                        "effect_allele_FLIP",
                        "other_allele_FLIP",
                    ]
                ).cast(str),
                pl.lit(True).alias("ambiguous"),
            ]
        ).pipe(complement_valid_alleles, ["REF"])
    ).with_columns(
        pl.when(pl.col("REF_FLIP") == pl.col("ALT"))
        .then(pl.col("ambiguous"))
        .otherwise(False)
    )

    if remove_ambiguous:
        logger.debug("Labelling ambiguous variants with exclude flag")
        return (
            ambig.with_columns(
                pl.when(pl.col("ambiguous"))
                .then(True)
                .otherwise(False)
                .alias("exclude_ambiguous")
            )
            .with_columns(max=pl.max_horizontal("exclude", "exclude_ambiguous"))
            .drop(["exclude", "exclude_ambiguous"])
            .rename({"max": "exclude"})
        )
    else:
        return (
            ambig.with_column(pl.lit(False).alias("exclude_ambiguous"))
            .with_columns(max=pl.max_horizontal("exclude", "exclude_ambiguous"))
            .drop(["exclude", "exclude_ambiguous"])
            .rename({"max": "exclude"})
        )


def _label_multiallelic(df: pl.LazyFrame, remove_multiallelic: bool) -> pl.LazyFrame:
    """Label multiallelic variants with exclude flag

    (Multiallelic variants are already labelled with the "is_multiallelic" column in match.preprocess)
    """
    if remove_multiallelic:
        logger.debug("Labelling multiallelic matches with exclude flag")
        return df.with_columns(
            pl.when(pl.col("is_multiallelic"))
            .then(True)
            .otherwise(pl.col("exclude"))  # don't overwrite existing exclude flags
            .alias("exclude")
        )
    else:
        logger.debug("Not excluding multiallelic variants")
        return df


def _label_flips(df: pl.LazyFrame, skip_flip: bool) -> pl.LazyFrame:
    df = df.with_columns(
        pl.when(pl.col("match_type").str.contains("_flip"))
        .then(True)
        .otherwise(False)
        .alias("match_flipped")
    )
    if skip_flip:
        logger.debug("Labelling flipped matches with exclude flag")
        return df.with_columns(
            pl.when(pl.col("match_flipped"))
            .then(True)
            .otherwise(pl.col("exclude"))  # don't overwrite existing exclude flags
            .alias("exclude")
        )
    else:
        logger.debug("Not excluding flipped matches")
        return df.with_columns(match_IDs=pl.lit("NA"))


def _label_filter(df: pl.LazyFrame, filter_IDs: list) -> pl.LazyFrame:
    if filter_IDs is not None:
        nIDs = len(filter_IDs)
        logger.debug(
            "Excluding variants that are not in ID list (read {} IDs)".format(nIDs)
        )
        df = df.with_columns(pl.col("ID").is_in(filter_IDs).alias("match_IDs"))
        return df.with_columns(
            pl.when(pl.col("match_IDs"))
            .then(True)
            .otherwise(pl.col("exclude"))
            .alias("exclude")
        )
    else:
        return df
