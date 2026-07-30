"""Target genome shapes and CLI drivers shared by the integration tests.

A shard varies along two axes at once (ADR 0003): the variant axis (which
chromosome) and the sample axis (which sample batch). Every fixture here is one
point in that space, built from the same records so they all hold the same
(sample, variant) pairs and must therefore produce the same scores:

- ``unsplit_target``          one file, no split at all (the ALL sentinel)
- ``chromosome_targets``      variant axis only, one file per chromosome
- ``sample_batch_targets``    sample axis only, two whole-genome batches
- ``grid_targets``            both axes, chromosome x sample batch
- ``variant_split_targets``   two files per chromosome sharing every sample, which
                              docs/tickets/0005 rejects
"""

import gzip
import pathlib
import sys

import polars as pl
import pysam
import pytest

from pgscatalog.calc.cli.main import main

# two chromosomes, so the variant axis is exercised without a slow fixture
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


def _subset_samples(
    source: pathlib.Path, path: pathlib.Path, samples: list[str]
) -> pathlib.Path:
    """Slice a VCF by sample, which is what a sample batch is"""
    pysam.bcftools.view(
        "-s",
        ",".join(samples),
        "-Oz",
        "-o",
        str(path),
        str(source),
        catch_stdout=False,
    )
    pysam.bcftools.index(str(path))
    return path


@pytest.fixture(scope="session")
def target_records(vcf_path):
    """Simple biallelic SNPs from two chromosomes, plus the header to write them"""
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

    for chrom, records in by_chrom.items():
        assert len(records) > 10, f"not enough test variants on {chrom}"

    return header, by_chrom


@pytest.fixture(scope="session")
def target_samples(target_records) -> list[str]:
    header, _ = target_records
    return list(header.samples)


@pytest.fixture(scope="session")
def scorefile(target_records, tmp_path_factory) -> pathlib.Path:
    """A scoring file covering every variant in the target genomes"""
    _, by_chrom = target_records
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


@pytest.fixture(scope="session")
def n_score_variants(target_records) -> int:
    _, by_chrom = target_records
    return sum(len(records) for records in by_chrom.values())


@pytest.fixture(scope="session")
def unsplit_target(target_records, tmp_path_factory) -> pathlib.Path:
    """Whole-genome input, which lands under the ALL sentinel chromosome"""
    header, by_chrom = target_records
    path = tmp_path_factory.mktemp("unsplit") / "unsplit.vcf.gz"
    _write_vcf(path, header, [r for recs in by_chrom.values() for r in recs])
    return path


@pytest.fixture(scope="session")
def chromosome_targets(target_records, tmp_path_factory) -> list[tuple]:
    """One file per chromosome: shards on the variant axis only"""
    header, by_chrom = target_records
    out_dir = tmp_path_factory.mktemp("by_chromosome")

    targets = []
    for chrom, records in by_chrom.items():
        path = out_dir / f"chr{chrom}.vcf.gz"
        _write_vcf(path, header, records)
        targets.append((path, chrom))

    return targets


@pytest.fixture(scope="session")
def sample_batch_targets(unsplit_target, target_samples, tmp_path_factory) -> list:
    """Two whole-genome files holding disjoint samples: the sample axis only"""
    out_dir = tmp_path_factory.mktemp("sample_batches")
    return [
        _subset_samples(unsplit_target, out_dir / f"batch{i}.vcf.gz", half)
        for i, half in enumerate((target_samples[::2], target_samples[1::2]), start=1)
    ]


@pytest.fixture(scope="session")
def grid_targets(chromosome_targets, target_samples, tmp_path_factory) -> list[tuple]:
    """Chromosome x sample batch, the shape docs/tickets/0008 has to reduce"""
    out_dir = tmp_path_factory.mktemp("grid")

    targets = []
    for source, chrom in chromosome_targets:
        for i, half in enumerate(
            (target_samples[::2], target_samples[1::2]), start=1
        ):
            path = _subset_samples(
                source, out_dir / f"chr{chrom}.batch{i}.vcf.gz", half
            )
            targets.append((path, chrom))

    return targets


@pytest.fixture(scope="session")
def variant_split_targets(target_records, tmp_path_factory) -> list[tuple]:
    """Two files per chromosome sharing every sample, splitting the variants.

    docs/tickets/0005 makes this an error: leaves under one chromosome must have
    disjoint samples, because the chromosome *is* the variant axis.
    """
    header, by_chrom = target_records
    out_dir = tmp_path_factory.mktemp("variant_split")

    targets = []
    for chrom, records in by_chrom.items():
        for i, half in enumerate((records[::2], records[1::2]), start=1):
            path = out_dir / f"chr{chrom}.half{i}.vcf.gz"
            _write_vcf(path, header, half)
            targets.append((path, chrom))

    return targets


@pytest.fixture(scope="session")
def golden_scores() -> pl.DataFrame:
    """Scores for the fixture data, generated before SCORE was split into stages.

    See tests/data/golden/README.md for how to regenerate.
    """
    path = pathlib.Path(__file__).parent.parent / "data" / "golden" / "scores.txt.gz"
    return pl.read_csv(path)


def run_cli(argv: list[str], monkeypatch) -> None:
    """Drive the CLI exactly as a Nextflow task does"""
    monkeypatch.setattr(sys, "argv", argv)
    main()


def load_target(
    target: str, cache_dir: pathlib.Path, scorefile, monkeypatch
) -> pathlib.Path:
    """LOAD one target genome into its own cache, like a Nextflow task would

    ``target`` is PATH or PATH:CHROM.
    """
    run_cli(
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


def load_targets(
    targets: list, cache_dir: pathlib.Path, scorefile, monkeypatch
) -> list[pathlib.Path]:
    """LOAD every shard into its own cache, one task each"""
    zips = []
    for i, target in enumerate(targets):
        if isinstance(target, tuple):
            path, chrom = target
            spec = f"{path}:{chrom}"
        else:
            spec = str(target)
        zips.append(
            load_target(spec, cache_dir / f"cache{i}", scorefile, monkeypatch)
        )
    return zips


def run_pipeline(
    targets: list,
    work_dir: pathlib.Path,
    scorefile,
    monkeypatch,
    min_overlap: str = "0.75",
    score_paths: list | None = None,
    max_memory_gb: str = "16",
    threads: str = "1",
) -> pl.DataFrame:
    """Drive MATCH -> MATCH_COMBINE -> SCORE -> SCORE_COMBINE, mirroring the
    Nextflow DAG: one task per chromosome, then one, then one per shard, then one.

    Every stage is its own task with its own resource limits (ADR 0001), so the
    limits are passed to each one just like Nextflow does.
    """
    score_paths = [str(p) for p in (score_paths or [scorefile])]
    resources = ["--max_memory_gb", max_memory_gb, "--threads", threads]
    zips = load_targets(targets, work_dir / "load", scorefile, monkeypatch)

    # MATCH fans out per chromosome, so group the shards by the chromosome they
    # were loaded under
    by_chromosome: dict[str, list[pathlib.Path]] = {}
    for target, zip_path in zip(targets, zips, strict=True):
        chrom = target[1] if isinstance(target, tuple) else "ALL"
        by_chromosome.setdefault(chrom, []).append(zip_path)

    matches, samples = [], []
    for chrom, chrom_zips in by_chromosome.items():
        out_dir = work_dir / "match" / chrom
        run_cli(
            [
                "pgsc_calc",
                "match",
                "--zarr_zip_file",
                *[str(z) for z in chrom_zips],
                "--score_paths",
                *score_paths,
                "--out_dir",
                str(out_dir),
                *resources,
            ],
            monkeypatch,
        )
        matches.append(out_dir / "matches.parquet")
        samples.append(out_dir / "samples.parquet")

    combine_dir = work_dir / "match_combine"
    run_cli(
        [
            "pgsc_calc",
            "match-combine",
            "--matches",
            *[str(p) for p in matches],
            "--samples",
            *[str(p) for p in samples],
            "--score_paths",
            *score_paths,
            "--sampleset",
            SAMPLESET,
            "--min_overlap",
            min_overlap,
            "--out_dir",
            str(combine_dir),
            *resources,
        ],
        monkeypatch,
    )
    match_table = combine_dir / "match_table.parquet"

    partials = []
    for i, zip_path in enumerate(zips):
        out_dir = work_dir / "score" / str(i)
        run_cli(
            [
                "pgsc_calc",
                "score",
                "--zarr_zip_file",
                str(zip_path),
                "--score_paths",
                *score_paths,
                "--match_table",
                str(match_table),
                "--out_dir",
                str(out_dir),
                *resources,
            ],
            monkeypatch,
        )
        partials.append(out_dir / "partial_scores.parquet")

    final_dir = work_dir / "combine"
    run_cli(
        [
            "pgsc_calc",
            "combine",
            "--partial_scores",
            *[str(p) for p in partials],
            "--out_dir",
            str(final_dir),
            *resources,
        ],
        monkeypatch,
    )

    return pl.read_csv(final_dir / "scores.txt.gz")
