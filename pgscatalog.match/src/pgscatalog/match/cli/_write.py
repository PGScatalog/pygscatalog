""" This module contains internal functions for writing match results to scoring files and logs

It expects Config class attributes to be set up before being called
"""
import gzip
import itertools

from ..lib.plinkscorefiles import PlinkScoreFiles

from ._config import Config


def write_matches(matchresults, score_df):
    """Write matchresults out to scoring files and logs"""
    match (Config.SPLIT, Config.COMBINED):
        case (True, True):
            # requires extra work: first write split
            outfs = matchresults.write_scorefiles(
                directory=Config.OUTDIR,
                split=True,
                score_df=score_df,
                min_overlap=Config.MIN_OVERLAP,
                **Config.MATCH_PARAMS,
            )
            # now re-combine without recomputing matches
            PlinkScoreFiles(*list(itertools.chain(*outfs))).merge(Config.OUTDIR)
        case (True, False) | (False, True):
            # split parameter can handle this case OK
            _ = matchresults.write_scorefiles(
                directory=Config.OUTDIR,
                split=Config.SPLIT,
                score_df=score_df,
                min_overlap=Config.MIN_OVERLAP,
                **Config.MATCH_PARAMS,
            )
        case _:
            raise ValueError

    write_log(matchresults=matchresults, score_df=score_df)
    # returns labelled and filtered data for checking after merging
    return matchresults.df


def write_log(matchresults, score_df):
    logfname = Config.OUTDIR / f"{Config.DATASET}_log.csv.gz"

    # summary log is smol
    matchresults.summary_log.write_csv(Config.OUTDIR / f"{Config.DATASET}_summary.csv")

    # this one can get big. gzip is slow, but everywhere
    with gzip.open(logfname, "wb") as f:
        matchresults.full_variant_log(score_df=score_df).collect().write_csv(f)
