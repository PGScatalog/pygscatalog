""" This module contains functions to filter match candidates and report overall match rates"""
import logging

import polars as pl

logger = logging.getLogger(__name__)


def filter_scores(
    scorefile: pl.LazyFrame, matches: pl.LazyFrame, min_overlap: float, dataset: str
) -> tuple[pl.LazyFrame, pl.DataFrame]:
    """Check overlap between filtered matches and scorefile

    Return a df with bad scores removed and a summary table"""
    filtered_matches: pl.LazyFrame = _filter_matches(matches)
    match_log: pl.LazyFrame = _join_filtered_matches(
        filtered_matches, scorefile, dataset
    ).with_columns(pl.col("best_match").fill_null(False))

    fail_rates: pl.DataFrame = _calculate_match_rate(
        match_log
    ).collect()  # collect for iteration

    scores: list[pl.DataFrame] = []
    for accession, rate in zip(
        fail_rates["accession"].to_list(), fail_rates["fail_rate"].to_list()
    ):
        if rate <= (1 - min_overlap):
            df: pl.DataFrame = pl.DataFrame(
                {
                    "accession": [accession],
                    "score_pass": [True],
                    "match_rate": [1 - rate],
                }
            )
            logger.debug(
                f"Score {accession} passes minimum matching threshold ({1 - rate:.2%}  variants match)"
            )
            scores.append(df.with_columns(pl.col("accession").cast(pl.Categorical)))
        else:
            df: pl.DataFrame = pl.DataFrame(
                {
                    "accession": [accession],
                    "score_pass": [False],
                    "match_rate": [1 - rate],
                }
            )
            logger.error(
                f"Score {accession} fails minimum matching threshold ({1 - rate:.2%} variants match)"
            )
            scores.append(df.with_columns(pl.col("accession").cast(pl.Categorical)))

    score_summary: pl.LazyFrame = pl.concat(scores).lazy()
    filtered_scores: pl.LazyFrame = filtered_matches.join(
        score_summary, on="accession", how="left"
    ).filter(pl.col("score_pass"))

    return filtered_scores, score_summary.collect()


def _calculate_match_rate(df: pl.LazyFrame) -> pl.LazyFrame:
    logger.debug("Calculating overlap between target genome and scoring file")
    return (
        df.group_by("accession")
        .agg([pl.all().len(), (pl.col("match_type").is_null()).sum().alias("no_match")])
        .with_columns((pl.col("no_match") / pl.col("count")).alias("fail_rate"))
    )


def _filter_matches(df: pl.LazyFrame) -> pl.LazyFrame:
    logger.debug("Filtering to best_match variants (with exclude flag = False)")
    return df.filter((pl.col("best_match")) & (~pl.col("exclude")))


def _join_filtered_matches(
    matches: pl.LazyFrame, scorefile: pl.LazyFrame, dataset: str
) -> pl.LazyFrame:
    return (
        scorefile.join(matches, on=["row_nr", "accession"], how="left")
        .with_columns(pl.lit(dataset).cast(pl.Categorical).alias("dataset"))
        .select(pl.exclude("^.*_right$"))
    )


def _match_keys() -> list[str]:
    return [
        "chr_name",
        "chr_position",
        "effect_allele",
        "other_allele",
        "accession",
        "effect_type",
        "effect_weight",
    ]
