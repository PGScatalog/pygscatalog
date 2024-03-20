import csv
import functools
import gzip
import itertools
import logging
import operator
import pathlib
from collections import namedtuple

from xopen import xopen

logger = logging.getLogger(__name__)


def _read_map(in_map, col_from, col_to):
    """Read a column from a mapping file into a dictionary"""
    mapping = {}
    logger.debug(f"Reading map {in_map}")
    with xopen(in_map) as map_f:
        reader = csv.DictReader(map_f, delimiter=" ")
        for line in reader:
            mapping[line[col_from]] = line[col_to]

    return mapping


def _read_maps(map_list, col_from, col_to):
    """Read mapping files and return one big dictionary"""
    maps = (_read_map(in_map=x, col_from=col_from, col_to=col_to) for x in map_list)
    return functools.reduce(operator.ior, maps)


def _get_outf_path(chrom, dataset):
    return pathlib.Path(f"{dataset}_{chrom}_relabelled.gz")


def _grab_header(in_target, comment_char):
    """Grab header from a TSV file which contains comment lines that are skipped"""
    for line in in_target:
        if line.startswith(comment_char):
            continue
        return line.strip().split()


def _relabel(in_path, mapping, target_col, comment_char):
    """Revalue a target column in a file based on the mapping dictionary"""
    with xopen(in_path) as in_f:
        h = _grab_header(in_f, comment_char)
        reader = csv.DictReader(in_f, delimiter="\t", fieldnames=h)
        for line in reader:
            line[target_col] = mapping[line[target_col]]
            yield line


def _parse_chrom(x):
    """Grab chromosome from the ID key of a dictionary

    (Some file types may not contain a CHROM column)"""
    return x["ID"].split(":")[0]


def relabel_write(relabelled, dataset, split_output, combined_output, out_dir):
    """Write relabelled data to a new file, optionally splitting by chromosome

    >>> from ._config import Config
    >>> import tempfile, os
    >>> maps = [Config.ROOT_DIR / "tests" / "data"/ "relabel_map_chr1.txt.gz", Config.ROOT_DIR / "tests" / "data" / "relabel_map_chr8.txt.gz"]
    >>> in_target = Config.ROOT_DIR / "tests" / "data" / "hgdp_ALL_additive_0.scorefile"
    >>> args = RelabelArgs()
    >>> x = relabel(in_path=in_target, map_paths=maps, relabel_args=args)
    >>> with tempfile.TemporaryDirectory() as tmp_dir:
    ...     relabel_write(x, dataset="test", split_output=True, combined_output=True, out_dir=tmp_dir)
    ...     sorted(os.listdir(tmp_dir))
    ['test_1_relabelled.gz', 'test_8_relabelled.gz', 'test_ALL_relabelled.gz']
    """
    if combined_output:
        combined_output = out_dir / pathlib.Path(
            _get_outf_path(chrom="ALL", dataset=dataset)
        )
        if combined_output.exists():
            raise FileExistsError(f"{combined_output}")

    for i, (chrom, group) in enumerate(itertools.groupby(relabelled, _parse_chrom)):
        # grab a batch of relabelled lines, grouped by chromosome
        batch = list(group)
        fieldnames = batch[0].keys()

        if split_output:
            outf_path = out_dir / _get_outf_path(chrom=chrom, dataset=dataset)
            logger.debug(f"Writing chrom {chrom} to {outf_path}")

            with gzip.open(outf_path, "wt") as csv_file:
                writer = csv.DictWriter(csv_file, delimiter="\t", fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(batch)

        if combined_output:
            logger.debug(f"Appending to {combined_output}")
            with gzip.open(combined_output, "at") as csv_file:
                writer = csv.DictWriter(csv_file, delimiter="\t", fieldnames=fieldnames)
                if i == 0:
                    writer.writeheader()

                writer.writerows(batch)


RelabelArgs = namedtuple(
    "RelabelArgs",
    ["comment_char", "dataset", "target_col", "map_col_from", "map_col_to"],
    defaults=["##", "DATASET", "ID", "ID_TARGET", "ID_REF"],
)


def relabel(*, in_path, map_paths, relabel_args):
    """Revalue a target column in a file based on mapping files

    Returns a generator of relabelled rows (one row = one dictionary).

    >>> from ._config import Config
    >>> relabel_args = RelabelArgs()
    >>> maps = [Config.ROOT_DIR / "tests" / "data" / "relabel_map_chr1.txt.gz", Config.ROOT_DIR / "tests" / "data" / "relabel_map_chr8.txt.gz"]

    The maps show that the variant on chromosome 8 has different ID across target and reference:

    >>> _read_maps(maps, col_from=relabel_args.map_col_from, col_to=relabel_args.map_col_to)
    {'1:11796321:G:A': '1:11796321:G:A', '8:127401060:G:T': '8:127401060:T:G'}

    >>> in_target = Config.ROOT_DIR / "tests" / "data" / "hgdp_ALL_additive_0.scorefile"
    >>> with open(in_target) as f:
    ...     for x in f:
    ...         x.split()
    ['ID', 'effect_allele', 'PGS000802_hmPOS_GRCh38']
    ['1:11796321:G:A', 'A', '0.16']
    ['8:127401060:G:T', 'G', '0.217']

    The variant on chromosome 8 will get relabelled:

    >>> x = relabel(in_path=in_target, map_paths=maps, relabel_args=relabel_args)
    >>> for row in x:
    ...     row
    {'ID': '1:11796321:G:A', 'effect_allele': 'A', 'PGS000802_hmPOS_GRCh38': '0.16'}
    {'ID': '8:127401060:T:G', 'effect_allele': 'G', 'PGS000802_hmPOS_GRCh38': '0.217'}
    """
    mapping = _read_maps(
        map_list=map_paths,
        col_from=relabel_args.map_col_from,
        col_to=relabel_args.map_col_to,
    )
    yield from _relabel(
        in_path=in_path,
        mapping=mapping,
        target_col=relabel_args.target_col,
        comment_char=relabel_args.comment_char,
    )
