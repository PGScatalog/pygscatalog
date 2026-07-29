"""Acceptance tests for the nested /{sampleset}/{chromosome}/{filename} cache
hierarchy (docs/tickets/0001, ADR 0003).

The interesting shapes are:

- more than one file per chromosome, which the old depth-2 layout couldn't store
- unsplit whole-genome input, which lands under the ALL sentinel chromosome

Both are scored end-to-end and must agree, because they hold the same
(sample, variant) pairs.
"""

import gzip
import pathlib
import sys
import zipfile

import polars as pl
import pysam
import pytest
import zarr
import zarr.storage

from pgscatalog.calc.cli.load import unzip_zarr
from pgscatalog.calc.cli.main import main

# two chromosomes, split into two files each
TEST_CHROMS = ("1", "8")
SAMPLESET = "test"


def _is_ambiguous(ref: str, alt: str) -> bool:
    """Ambiguous variants get excluded by the matcher, which would drag the match
    rate down for no reason here"""
    return {ref, alt} in ({"A", "T"}, {"C", "G"})


def _write_vcf(path: pathlib.Path, header: pysam.VariantHeader, records: list) -> None:
    with pysam.VariantFile(str(path), mode="w", header=header) as out:
        for rec in records:
            out.write(rec)
    pysam.bcftools.index(str(path))


@pytest.fixture(scope="module")
def targets(vcf_path, tmp_path_factory):
    """Simple biallelic SNPs from two chromosomes, written out twice:

    - as four files (2 chromosomes x 2 files), the nested case
    - as one whole-genome file, the ALL sentinel case
    """
    out_dir = tmp_path_factory.mktemp("targets")
    by_chrom: dict[str, list] = {chrom: [] for chrom in TEST_CHROMS}

    with pysam.VariantFile(str(vcf_path)) as vcf:
        header = vcf.header.copy()
        for rec in vcf.fetch():
            if rec.chrom not in by_chrom:
                continue
            if len(rec.ref) != 1 or len(rec.alts) != 1 or len(rec.alts[0]) != 1:
                continue
            if _is_ambiguous(rec.ref, rec.alts[0]):
                continue
            by_chrom[rec.chrom].append(rec)

    # 2 files per chromosome, splitting the variants between them
    split: list[tuple[pathlib.Path, str]] = []
    for chrom, records in by_chrom.items():
        assert len(records) > 10, f"not enough test variants on {chrom}"
        for batch, half in enumerate((records[::2], records[1::2]), start=1):
            path = out_dir / f"chr{chrom}.batch{batch}.vcf.gz"
            _write_vcf(path, header, half)
            split.append((path, chrom))

    unsplit = out_dir / "unsplit.vcf.gz"
    _write_vcf(unsplit, header, [r for recs in by_chrom.values() for r in recs])

    return split, unsplit, by_chrom


@pytest.fixture(scope="module")
def scorefile(targets, tmp_path_factory):
    """A scoring file covering every variant in the target genomes"""
    _, _, by_chrom = targets
    path = tmp_path_factory.mktemp("scorefile") / "PGS999998.txt.gz"

    with gzip.open(path, "wt") as f:
        f.write(
            "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\t"
            "effect_type\tis_duplicated\taccession\trow_nr\n"
        )
        row_nr = 0
        for records in by_chrom.values():
            for rec in records:
                weight = round(0.001 * (row_nr % 20) - 0.01, 6)
                f.write(
                    f"{rec.chrom}\t{rec.pos}\t{rec.alts[0]}\t{rec.ref}\t{weight}\t"
                    f"additive\tFalse\tPGS999998\t{row_nr}\n"
                )
                row_nr += 1

    return path


def _run(argv: list[str], monkeypatch) -> None:
    monkeypatch.setattr(sys, "argv", argv)
    main()


def _load(target: str, cache_dir: pathlib.Path, scorefile, monkeypatch) -> pathlib.Path:
    """LOAD one target genome into its own cache, like a Nextflow task would"""
    _run(
        [
            "pgsc_calc",
            "load",
            "--target_genomes",
            target,
            "--sampleset",
            SAMPLESET,
            "--format",
            "vcf",
            "--cache_dir",
            str(cache_dir),
            "--score_paths",
            str(scorefile),
        ],
        monkeypatch,
    )
    return cache_dir / "genotypes.zarr.zip"


def _score(
    zips: list[pathlib.Path], out_dir: pathlib.Path, scorefile, monkeypatch
) -> pl.DataFrame:
    _run(
        [
            "pgsc_calc",
            "score",
            "--zarr_zip_file",
            *[str(z) for z in zips],
            "--score_paths",
            str(scorefile),
            "--out_dir",
            str(out_dir),
            "--min_overlap",
            "0.75",
        ],
        monkeypatch,
    )
    return pl.read_csv(out_dir / "scores.txt.gz")


def test_two_files_per_chromosome(targets, scorefile, tmp_path, monkeypatch):
    """2 chromosomes x 2 files per chromosome loads, zips, and scores"""
    split, _, by_chrom = targets
    n_variants = sum(len(recs) for recs in by_chrom.values())

    zips = []
    for i, (path, chrom) in enumerate(split):
        zips.append(
            _load(f"{path}:{chrom}", tmp_path / f"cache{i}", scorefile, monkeypatch)
        )

    # each shard is stored at its own chromosome, so basenames no longer collide
    for (path, chrom), zip_path in zip(split, zips, strict=True):
        names = zipfile.ZipFile(zip_path).namelist()
        prefix = f"genotypes.zarr/{SAMPLESET}/{chrom}/{path.name}/"
        assert any(n.startswith(prefix + "genotypes/") for n in names), names[:5]

    scores = _score(zips, tmp_path / "out", scorefile, monkeypatch)

    assert scores["sampleset"].unique().to_list() == [SAMPLESET]
    assert scores["accession"].unique().to_list() == ["PGS999998"]
    # every sample gets one score
    assert scores.shape[0] == 3202
    # all four shards contributed, not just one
    assert scores["n_matched"].unique().to_list() == [n_variants]
    assert scores["score"].n_unique() > 1


def test_unsplit_input_uses_all_sentinel(targets, scorefile, tmp_path, monkeypatch):
    """Whole-genome input lands under ALL and scores the same as the split input"""
    split, unsplit, _ = targets

    unsplit_zip = _load(
        str(unsplit), tmp_path / "unsplit_cache", scorefile, monkeypatch
    )
    names = zipfile.ZipFile(unsplit_zip).namelist()
    prefix = f"genotypes.zarr/{SAMPLESET}/ALL/{unsplit.name}/"
    assert any(n.startswith(prefix + "genotypes/") for n in names), names[:5]

    unsplit_scores = _score(
        [unsplit_zip], tmp_path / "unsplit_out", scorefile, monkeypatch
    )

    split_zips = [
        _load(f"{path}:{chrom}", tmp_path / f"split{i}", scorefile, monkeypatch)
        for i, (path, chrom) in enumerate(split)
    ]
    split_scores = _score(split_zips, tmp_path / "split_out", scorefile, monkeypatch)

    # same (sample, variant) pairs however the cache is shaped, so same scores
    key = ["sampleset", "accession", "sample_id"]
    assert unsplit_scores.sort(key)["score"].to_list() == pytest.approx(
        split_scores.sort(key)["score"].to_list()
    )


def test_partial_extract_of_one_leaf(targets, scorefile, tmp_path, monkeypatch):
    """Extracting a single leaf brings every ancestor zarr.json with it, whatever
    the depth of the hierarchy"""
    split, _, _ = targets
    path, chrom = split[0]
    zip_path = _load(f"{path}:{chrom}", tmp_path / "cache", scorefile, monkeypatch)

    extract_dir = tmp_path / "extracted"
    extract_dir.mkdir()
    group_path = f"genotypes.zarr/{SAMPLESET}/{chrom}/{path.name}"
    unzip_zarr(zip_path=zip_path, cache_path=extract_dir, group_path=group_path)

    # the hierarchy only parses if root, sampleset and chromosome metadata is present
    store = zarr.storage.LocalStore(root=extract_dir / "genotypes.zarr")
    root = zarr.open_group(store=store, mode="r")
    leaf = root[f"{SAMPLESET}/{chrom}/{path.name}"]
    assert "genotypes" in leaf

    # and nothing else got dragged along
    assert list(root[SAMPLESET].group_keys()) == [chrom]
