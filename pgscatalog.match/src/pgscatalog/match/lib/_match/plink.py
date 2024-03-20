import logging

import polars as pl

from pgscatalog.core import EffectType

logger = logging.getLogger(__name__)


def plinkify(df):
    """Prepare dataframe to align with plink2 expectations for a --scorefile

    Split effect types: weight calculation differs by effect type
    Variants with the same ID and different effect alleles must be split into different files

    Input df must contain match candidates filtered to best matched variant
    """
    min_cols = [
        "accession",
        "effect_type",
        "chr_name",
        "ID",
        "matched_effect_allele",
        "effect_weight",
    ]
    df = df.select(min_cols).lazy()

    # 1. split by effect type
    effect_frames = zip(
        [EffectType.ADDITIVE, EffectType.DOMINANT, EffectType.RECESSIVE],
        _split_effect_type(df),
        strict=True,
    )

    # 2. deduplicate
    deduped = {}
    for effect_type, effect_frame in effect_frames:
        deduped.setdefault(effect_type, [])
        deduped[effect_type].append(effect_frame)

    return deduped


def pivot_score(df: pl.DataFrame) -> pl.DataFrame:
    """Format a dataframe to plink2 --score standard
    Minimum example:
    ID | effect_allele | effect_weight
    Multiple scores are OK too:
    ID | effect_allele | weight_1 | ... | weight_n
    """
    return (
        df.pivot(
            index=["ID", "matched_effect_allele", "effect_type"],
            values="effect_weight",
            columns="accession",
        )
        .rename({"matched_effect_allele": "effect_allele"})
        .fill_null(strategy="zero")
        .drop("effect_type")
    )


def _check_column_types(matches: pl.LazyFrame):
    logger.debug("Checking column types")
    # these columns are most important for writing out
    correct_schema = {
        "chr_name": pl.Categorical,
        "chr_position": pl.UInt64,
        "ID": pl.Utf8,
        "matched_effect_allele": pl.Categorical,
        "effect_weight": pl.Float64,
        "effect_type": pl.Categorical,
        "accession": pl.Categorical,
    }
    col_types = {
        x: matches.schema.get(x)
        for x in list((matches.schema.keys() & correct_schema.keys()))
    }
    if not (col_types == correct_schema):
        logger.critical(
            f"MISMATCHED SCHEMA\nCurrent columns: {col_types}\nCorrect schema:{correct_schema}"
        )
        raise ValueError


def _split_effect_type(
    df: pl.LazyFrame,
) -> tuple[pl.LazyFrame, pl.LazyFrame, pl.LazyFrame]:
    additive = df.filter(pl.col("effect_type") == "additive")
    dominant = df.filter(pl.col("effect_type") == "dominant")
    recessive = df.filter(pl.col("effect_type") == "recessive")
    return additive, dominant, recessive


def _deduplicate_variants(df: pl.LazyFrame) -> list[pl.LazyFrame | None]:
    """Find variant matches that have duplicate identifiers
    When merging a lot of scoring files, sometimes a variant might be duplicated
    this can happen when the matched effect allele differs at the same position, e.g.:
        - chr1: chr2:20003:A:C A 0.3 NA
        - chr1: chr2:20003:A:C C NA 0.7
    where the last two columns represent different scores.  plink demands
    unique identifiers! so need to split, score, and sum later
    Parameters:
    df: A dataframe containing all matches, with columns ID, effect_allele, and
        effect_weight
    Returns:
        A list of dataframes, with unique ID - matched effect allele combinations
    """
    if df.select("ID").head().collect().is_empty():
        logger.info("Empty input: skipping deduplication")
        return [None]
    else:
        logger.debug("Deduplicating variants")

    # 1. unique ID - EA is important because normal duplicates are already
    #   handled by pivoting, and it's pointless to split them unnecessarily
    # 2. use cumcount to number duplicate IDs
    # 3. join cumcount data on original DF, use this data for splitting
    # note: effect_allele should be equivalent to matched_effect_allele
    ea_count: pl.LazyFrame = (
        df.select(["ID", "matched_effect_allele"])
        .unique()
        .with_columns(
            [
                pl.col("ID").cumcount().over(["ID"]).alias("cumcount"),
                pl.col("ID").count().over(["ID"]).alias("count"),
            ]
        )
    )

    dup_label: pl.LazyFrame = df.join(
        ea_count, on=["ID", "matched_effect_allele"], how="left"
    )

    # now split the matched variants, and make sure we don't lose any
    # cumcount = ngroup-1, so add 1
    n_splits: int = ea_count.select("cumcount").max().collect().item(0, 0) + 1
    df_lst: list = []

    for i in range(0, n_splits):
        x: pl.LazyFrame = dup_label.filter(pl.col("cumcount") == i).drop(
            ["cumcount", "count"]
        )
        if (df := x.collect()).is_empty():
            continue
        else:
            df_lst.append(df.lazy())

    return df_lst
