import logging

import polars as pl

from pgscatalog.core import EffectType

logger = logging.getLogger(__name__)


def plinkify(df):
    """Prepare dataframe to align with plink2 expectations for a --scorefile

    Split effect types: weight calculation differs by effect type
    Variants with the same ID and different effect alleles must be split into different files

    Input df must contain match candidates filtered to best matched variant

    Some fake match data where effect_allele is the same but other_allele is different:

    >>> data = {'row_nr': [1, 2, 3], 'chr_name': ['11', '11', '11'], 'chr_position': [69331418, 69331418, 69331419], 'effect_allele': ['A', 'C', 'T'], 'other_allele': ['C', 'C', 'T'], 'effect_weight': ['1', '1.5', '2'], 'effect_type': ['additive', 'additive', 'additive'], 'accession': ['pgs1', 'pgs2', 'pgs3'],'ID': ['11:69331418:A:C', '11:69331418:A:C', '11:69331419:T:T'],'matched_effect_allele': ['A', 'C', 'T']}
    >>> plinked = plinkify(pl.DataFrame(data).lazy())
    >>> dfs = sorted((x.collect() for x in plinked[EffectType.ADDITIVE]), key= lambda x: len(x))
    >>> len(dfs)
    2

    >>> dfs[0].select(["row_nr", "accession", "ID", "matched_effect_allele"])
    shape: (1, 4)
    ┌────────┬───────────┬─────────────────┬───────────────────────┐
    │ row_nr ┆ accession ┆ ID              ┆ matched_effect_allele │
    │ ---    ┆ ---       ┆ ---             ┆ ---                   │
    │ i64    ┆ str       ┆ str             ┆ str                   │
    ╞════════╪═══════════╪═════════════════╪═══════════════════════╡
    │ 2      ┆ pgs2      ┆ 11:69331418:A:C ┆ C                     │
    └────────┴───────────┴─────────────────┴───────────────────────┘

    >>> dfs[1].select(["row_nr", "accession", "ID", "matched_effect_allele"])
    shape: (2, 4)
    ┌────────┬───────────┬─────────────────┬───────────────────────┐
    │ row_nr ┆ accession ┆ ID              ┆ matched_effect_allele │
    │ ---    ┆ ---       ┆ ---             ┆ ---                   │
    │ i64    ┆ str       ┆ str             ┆ str                   │
    ╞════════╪═══════════╪═════════════════╪═══════════════════════╡
    │ 1      ┆ pgs1      ┆ 11:69331418:A:C ┆ A                     │
    │ 3      ┆ pgs3      ┆ 11:69331419:T:T ┆ T                     │
    └────────┴───────────┴─────────────────┴───────────────────────┘

    When merging a lot of scoring files, sometimes a variant might be duplicated
    this can happen when the matched effect allele differs at the same position, e.g.:
        - chr1: chr2:20003:A:C A 0.3 NA
        - chr1: chr2:20003:A:C C NA 0.7
    where the last two columns represent different scores.  plink demands
    unique identifiers! so need to split, score, and sum later.
    """
    min_cols = [
        "row_nr",
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
        deduped[effect_type] = _deduplicate_variants(effect_frame, effect_type)

    return deduped


def pivot_score(df: pl.DataFrame) -> pl.DataFrame:
    """Format a dataframe to plink2 --score standard
    Minimum example:
    ID | effect_allele | effect_weight
    Multiple scores are OK too:
    ID | effect_allele | weight_1 | ... | weight_n
    """
    # fill_null: value = 0 not strategy = "zero"
    # we always handle effect weights as strings because we never modify them
    return (
        df.pivot(
            index=["ID", "matched_effect_allele", "effect_type"],
            values="effect_weight",
            columns="accession",
        )
        .rename({"matched_effect_allele": "effect_allele"})
        .fill_null(value="0")
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


def _deduplicate_variants(df: pl.LazyFrame, effect_type: str) -> list[pl.LazyFrame]:
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
    if df.select(pl.len()).collect().item(0, 0) == 0:
        logger.info(f"No variants with {effect_type=}, skipping deduplication")
        return [df]

    logger.debug(
        f"Splitting {effect_type} variants with different effect alleles and the same ID"
    )

    # 0: make sure input is sorted for test output to be consistent
    df = df.sort("row_nr")

    # 1. count the number of unique matched alleles per ID
    id_alleles_count: pl.DataFrame = (
        df.unique(["matched_effect_allele", "ID"], maintain_order=True)
        .with_columns(
            pl.cum_count("matched_effect_allele").over("ID").alias("cum_count")
        )
        .collect()
    )

    # 2: use the cumulative count to create groups to iterate over when splitting
    groups = id_alleles_count.group_by(["cum_count"], maintain_order=True)

    # 3. unique() will have dropped rows. it's important to keep them
    dropped_variants: pl.LazyFrame = df.filter(
        ~pl.col("row_nr").is_in(id_alleles_count.select("row_nr"))
    )

    ldf_lst = []
    for i, (_, x) in enumerate(groups):
        tempdf: pl.LazyFrame = df.filter(
            pl.col("row_nr").is_in(x.select("row_nr").to_series())
        ).sort("row_nr")
        if i == 0:
            logger.info("Initialising deduplicated dataframe")
            tempdf: pl.LazyFrame = pl.concat([dropped_variants, tempdf]).sort("row_nr")
        elif i > 0:
            logger.info(f"Duplicate variants detected, creating dataframe {i}")
        ldf_lst.append(tempdf)

    # note: this function intentionally avoided joining, variants might go missing
    # let's double-check the number of variants at the beginning and end
    original_length: int = df.select(pl.len()).collect().item(0, 0)
    new_length: int = sum(x.select(pl.len()).collect().item(0, 0) for x in ldf_lst)

    if original_length != new_length:
        raise ValueError(
            f"Some variants went missing during deduplication ({original_length} != {new_length}) "
        )

    logger.info("Split dataframe lengths are consistent with input after deduplicating")
    return ldf_lst
