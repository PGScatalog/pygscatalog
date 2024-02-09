""" This module contains a simple class that stores CLI configuration as class attributes """


class Config:
    # what's the label for the target variants data?
    DATASET = None
    # where do scorefiles get written to?
    OUTDIR = None
    # where do we keep stable-ish arrow files? (matchresults)
    MATCHTMP = None
    # where do we keep unstable arrow files
    TMPDIR = None
    # filter variant matches to a single chromosome?
    CHROM = None
    # scorefile path
    SCOREFILE = None
    # clean up temporary arrow files created by context manager?
    CLEANUP = False
    # scorefile output options
    SPLIT = False
    COMBINED = False
    # matching parameter arguments
    MATCH_KWARGS = [
        "keep_first_match",
        "remove_ambiguous",
        "skip_flip",
        "filter_IDs",
        "remove_multiallelic",
    ]
    MATCH_PARAMS = None
    MIN_OVERLAP = 0.75
