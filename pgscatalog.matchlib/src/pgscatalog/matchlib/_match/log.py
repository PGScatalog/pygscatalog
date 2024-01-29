import logging

import polars as pl

logger = logging.getLogger(__name__)


def make_logs(scorefile: pl.LazyFrame, match_candidates: pl.LazyFrame, dataset: str):
    # summary log -> aggregated from best matches (one per scoring file line)
    # big log -> unaggregated, written to compressed gzip, possibly multiple matches per scoring file line
    big_log = _join_match_candidates(
        scorefile=scorefile, matches=match_candidates, dataset=dataset
    )

    return _prettify_log(big_log)


def make_summary_log(
    match_candidates: pl.LazyFrame,
    scorefile: pl.LazyFrame,
    filter_summary: pl.LazyFrame,
    dataset: str,
) -> pl.LazyFrame:
    """Make an aggregated table"""
    logger.debug("Aggregating best matches into a summary table")
    best_matches: pl.LazyFrame = match_candidates.filter(pl.col("best_match"))
    return (
        scorefile.join(best_matches, on=["row_nr", "accession"], how="outer")
        .select(pl.exclude("^.*_right$"))
        .with_columns(
            [
                pl.col("match_status").fill_null(value="unmatched"),
                pl.lit(dataset).alias("dataset"),
            ]
        )  # fill in unmatched variants
        .groupby(
            [
                "dataset",
                "accession",
                "ambiguous",
                "is_multiallelic",
                "match_flipped",
                "duplicate_best_match",
                "duplicate_ID",
                "match_IDs",
                "match_status",
            ]
        )
        .agg(pl.all().len())
        .join(filter_summary, how="left", on="accession")
        .pipe(_prettify_summary)
    )


def check_log_count(scorefile: pl.LazyFrame, summary_log: pl.LazyFrame) -> None:
    """Check aggregated counts vs original from scoring file"""
    summary_count: pl.DataFrame = (
        summary_log.groupby(pl.col("accession")).agg(pl.sum("count"))
    ).collect()
    log_count: pl.DataFrame = (
        scorefile.groupby("accession")
        .agg(pl.all().len())
        .collect()
        .join(summary_count, on="accession")
    )

    assert (
        log_count.get_column("count") == log_count.get_column("count_right")
    ).all(), "Log doesn't match input scoring file"
    logger.debug("Log matches input scoring file")


def _prettify_summary(df: pl.LazyFrame) -> pl.LazyFrame:
    keep_cols = [
        "dataset",
        "accession",
        "score_pass",
        "match_status",
        "ambiguous",
        "is_multiallelic",
        "duplicate_best_match",
        "duplicate_ID",
        "match_flipped",
        "match_IDs",
        "count",
        "percent",
    ]
    return df.with_columns(
        (pl.col("count") / pl.sum("count") * 100)
        .over(["dataset", "accession"])
        .alias("percent")
    ).select(keep_cols)


def _prettify_log(df: pl.LazyFrame) -> pl.LazyFrame:
    keep_cols = [
        "row_nr",
        "accession",
        "chr_name",
        "chr_position",
        "effect_allele",
        "other_allele",
        "effect_weight",
        "effect_type",
        "ID",
        "REF",
        "ALT",
        "matched_effect_allele",
        "match_type",
        "is_multiallelic",
        "ambiguous",
        "match_flipped",
        "best_match",
        "exclude",
        "duplicate_best_match",
        "duplicate_ID",
        "match_IDs",
        "match_status",
        "dataset",
    ]
    pretty_df = (
        df.select(keep_cols)
        .select(pl.exclude("^.*_right"))
        .sort(["accession", "row_nr", "chr_name", "chr_position", "match_status"])
    )
    return pretty_df


def _join_match_candidates(
    scorefile: pl.LazyFrame, matches: pl.LazyFrame, dataset: str
) -> pl.LazyFrame:
    """Join match candidates against the original scoring file"""
    logger.debug("Joining all match candidates against input scoring file")
    # make a raw log with all match candidates included
    raw_log = (
        scorefile.join(matches, on=["row_nr", "accession"], how="outer")
        .with_columns(pl.lit(dataset).alias("dataset").cast(pl.Categorical))
        .select(pl.exclude("^.*_right$"))
    ).with_columns(pl.col("match_status").fill_null("unmatched"))

    return raw_log
