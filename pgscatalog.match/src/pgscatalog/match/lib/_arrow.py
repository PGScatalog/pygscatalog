""" This module contains functions for streaming text files to arrow IPC files

It reimplements CSV reading using the polars API directly to be fast.

pgscatalog.core's implementation with csv.DictReader is slow with hundreds of millions of lines.
"""
import pathlib
import tempfile

import polars as pl
from functools import singledispatch

import polars.exceptions
from pgscatalog.core import NormalisedScoringFile, TargetVariants, TargetType


@singledispatch
def loose(path, tmpdir=None):
    """Loose an arrow :) Stream compressed text files into temporary arrow files

    This approach relies heavily on polars to read and write tabular data efficiently
    """
    # can currently only handle ScoringFile and TargetVariants
    raise NotImplementedError


@loose.register
def _(path: NormalisedScoringFile, tmpdir=None):
    """Write NormalisedScoringFiles to a list of arrow IPC files"""
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp()
    # alleles are first read as strings and converted to categoricals later
    dtypes = {
        "chr_name": pl.Categorical,
        "chr_position": pl.UInt64,
        "effect_allele": pl.String,
        "other_allele": pl.String,
        "effect_weight": pl.String,
        "effect_type": pl.Categorical,
        "is_duplicated": pl.Boolean,
        "accession": pl.Categorical,
        "row_nr": pl.UInt64,
    }

    reader = pl.read_csv_batched(path.path, dtypes=dtypes, separator="\t")
    return batch_read(reader, tmpdir=tmpdir, cols_keep=pl.all())


@loose.register
def _(path: TargetVariants, tmpdir=None):
    """Writes TargetVariants to a list of arrow IPC files"""
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp()

    match path.ftype:
        case TargetType.BIM:
            # no header
            dtypes = {
                "column_1": pl.Categorical,
                "column_2": pl.String,
                "column_3": pl.UInt8,
                "column_4": pl.UInt64,
                "column_5": pl.String,
                "column_6": pl.String,
            }
            new_colnames = ["#CHROM", "ID", "CM", "POS", "REF", "ALT"]
            header = False
            comment = None
        case TargetType.PVAR:
            # order of dict is important
            dtypes = {
                "#CHROM": pl.Categorical,
                "POS": pl.UInt64,
                "ID": pl.String,
                "REF": pl.String,
                "ALT": pl.String,
            }
            header = True
            comment = "##"
            new_colnames = None  # don't rename, use the header
        case _:
            raise NotImplementedError

    reader = pl.read_csv_batched(
        path.path,
        has_header=header,
        dtypes=dtypes,
        separator="\t",
        new_columns=new_colnames,
        comment_prefix=comment,
        batch_size=1000000,
    )
    cols_keep = ["#CHROM", "POS", "ID", "REF", "ALT"]
    return batch_read(reader, tmpdir=tmpdir, cols_keep=cols_keep)


def batch_read(reader, tmpdir, cols_keep) -> list[pathlib.Path]:
    """Read a CSV in batches and write them to temporary files"""
    arrowpaths = []
    # batch_size should be >= thread pool size, so tasks will be distributed amongst workers
    batch_size = 100
    batches = reader.next_batches(batch_size)
    while batches:
        arrowpath = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)
        df = pl.concat(batches).select(cols_keep)
        df.write_ipc(arrowpath)
        batches = reader.next_batches(batch_size)
        arrowpaths.append(pathlib.Path(arrowpath.name))

    return arrowpaths
