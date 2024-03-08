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
    """Make an aggregated table that contains the best match candidate for each row in the original scoring file"""
    logger.debug("Aggregating best matches into a summary table")
    cols = [
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

    best_matches: pl.LazyFrame = match_candidates.filter(pl.col("best_match"))

    return (
        scorefile.join(best_matches, on=["row_nr", "accession"], how="outer")
        .with_columns(
            pl.col("match_status").fill_null(value="unmatched"), dataset=pl.lit(dataset)
        )  # fill in unmatched variants
        .group_by(cols)
        .len()
        .rename({"len": "count"})
        .join(filter_summary, how="left", on="accession")
        .pipe(_prettify_summary)
    )


def check_log_count(scorefile: pl.LazyFrame, summary_log: pl.LazyFrame):
    """Check aggregated counts vs original from scoring file"""
    summary_count: pl.LazyFrame = summary_log.group_by("accession").agg(pl.sum("count"))

    log_count: pl.DataFrame = (
        scorefile.group_by("accession")
        .len()
        .rename({"len": "count"})
        .join(summary_count, on="accession")
        .collect()
    )

    for row in log_count.iter_rows():
        if (log_count := row[1]) != (score_count := row[2]):
            logger.critical("Variant log doesn't match input scoring file counts")
            raise ValueError(
                f"{row[0]} match failure {log_count=} doesn't match {score_count=}"
            )

    return True


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
    pretty_df = df.select(keep_cols).sort(
        ["accession", "row_nr", "chr_name", "chr_position", "match_status"]
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
