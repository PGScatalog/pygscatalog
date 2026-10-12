"""Cache shape invariants checked at MATCH and SCORE_COMBINE"""

import duckdb
import polars as pl
import pytest

from pgscatalog.calc.lib.cache.genotypecache import (
    Shard,
    check_absence_semantics,
    check_cache_shape,
    check_samples_disjoint,
    check_shards_are_rectangular,
)
from pgscatalog.calc.lib.score.match_combine import check_sample_union
from pgscatalog.calc.lib.score.score_combine import ScoreCombineStage


def _samples_parquet(tmp_path, name, rows):
    path = tmp_path / name
    pl.DataFrame(
        rows, schema=["sampleset", "chromosome", "sample_id"], orient="row"
    ).write_parquet(path)
    return path


def test_matching_sample_union_passes(tmp_path):
    paths = [
        _samples_parquet(
            tmp_path, f"chr{chrom}.parquet", [("test", chrom, s) for s in "abc"]
        )
        for chrom in ("1", "8")
    ]
    with duckdb.connect() as conn:
        check_sample_union(conn=conn, paths=paths)


def test_chromosome_missing_a_sample_is_rejected(tmp_path):
    paths = [
        _samples_parquet(tmp_path, "chr1.parquet", [("test", "1", s) for s in "abc"]),
        _samples_parquet(tmp_path, "chr8.parquet", [("test", "8", s) for s in "ab"]),
    ]
    with duckdb.connect() as conn, pytest.raises(ValueError, match=r"\['8'\]"):
        check_sample_union(conn=conn, paths=paths)


def test_sparse_scoring_file_does_not_trip_the_union_check(tmp_path):
    """A chromosome whose index query returned nothing is deleted from the cache
    and runs no MATCH task, so it simply isn't here to compare against"""
    paths = [
        _samples_parquet(tmp_path, "chr1.parquet", [("test", "1", s) for s in "abc"])
    ]
    with duckdb.connect() as conn:
        check_sample_union(conn=conn, paths=paths)


def test_shared_samples_are_rejected():
    shards = {
        "test/1/a.vcf.gz": {"samples": ["s1", "s2", "s3"]},
        "test/1/b.vcf.gz": {"samples": ["s3", "s4"]},
    }
    with pytest.raises(ValueError) as excinfo:
        check_samples_disjoint(shards)

    message = str(excinfo.value)
    assert "test/1/a.vcf.gz" in message
    assert "test/1/b.vcf.gz" in message
    assert "s3" in message


def test_disjoint_samples_pass():
    check_samples_disjoint(
        {
            "test/1/a.vcf.gz": {"samples": ["s1", "s2"]},
            "test/1/b.vcf.gz": {"samples": ["s3", "s4"]},
        }
    )


def test_unobserved_is_the_default_for_caches_without_a_declaration():
    check_absence_semantics({"test/1/a.vcf.gz": {"samples": []}})


def test_disagreeing_declarations_are_rejected():
    shards = {
        "test/1/a.vcf.gz": {"absence_semantics": "unobserved"},
        "test/1/b.vcf.gz": {"absence_semantics": "homref"},
    }
    with pytest.raises(ValueError, match="disagree about absence semantics"):
        check_absence_semantics(shards)


def test_homref_is_not_implemented():
    shards = {"test/1/a.vcf.gz": {"absence_semantics": "homref"}}
    with pytest.raises(NotImplementedError, match="not implemented yet"):
        check_absence_semantics(shards)


@pytest.fixture
def targetvariants():
    conn = duckdb.connect()
    conn.execute("""
        CREATE TABLE targetvariants (
            ref TEXT, variant_id TEXT, zarr_group TEXT
        );
    """)
    return conn


def test_rectangular_shards_pass(targetvariants):
    """Sentinel rows for queried-but-absent positions have a null ref, and are
    identical across shards by construction, so the check ignores them"""
    targetvariants.executemany(
        "INSERT INTO targetvariants VALUES (?, ?, ?)",
        [
            ("A", "1:100:A:C", "test/1/a.vcf.gz"),
            ("A", "1:200:A:G", "test/1/a.vcf.gz"),
            (None, "1:300::", "test/1/a.vcf.gz"),
            ("A", "1:100:A:C", "test/1/b.vcf.gz"),
            ("A", "1:200:A:G", "test/1/b.vcf.gz"),
            (None, "1:300::", "test/1/b.vcf.gz"),
        ],
    )
    check_shards_are_rectangular(
        targetvariants, {"test/1/a.vcf.gz": {}, "test/1/b.vcf.gz": {}}
    )


def test_ragged_shards_are_rejected_without_suggesting_homref(targetvariants):
    targetvariants.executemany(
        "INSERT INTO targetvariants VALUES (?, ?, ?)",
        [
            ("A", "1:100:A:C", "test/1/a.vcf.gz"),
            ("A", "1:200:A:G", "test/1/a.vcf.gz"),
            ("A", "1:100:A:C", "test/1/b.vcf.gz"),
        ],
    )
    with pytest.raises(ValueError) as excinfo:
        check_shards_are_rectangular(
            targetvariants, {"test/1/a.vcf.gz": {}, "test/1/b.vcf.gz": {}}
        )

    message = str(excinfo.value)
    assert "different variants" in message
    assert "(2 variants)" in message and "(1 variants)" in message
    # the offending variant is named, not just counted
    assert "1:200:A:G (only in test/1/a.vcf.gz)" in message
    # ragged shards are ambiguous in cause, so never advise declaring homref
    assert "homref" not in message.lower()


def test_shards_differing_by_swapped_variants_are_rejected():
    """Equal counts, different variants: comparing identifiers catches this, and
    comparing counts alone would not"""
    conn = duckdb.connect()
    conn.execute(
        "CREATE TABLE targetvariants (ref TEXT, variant_id TEXT, zarr_group TEXT);"
    )
    conn.executemany(
        "INSERT INTO targetvariants VALUES (?, ?, ?)",
        [
            ("A", "1:100:A:C", "test/1/a.vcf.gz"),
            ("A", "1:200:A:G", "test/1/a.vcf.gz"),
            ("A", "1:100:A:C", "test/1/b.vcf.gz"),
            ("A", "1:300:A:T", "test/1/b.vcf.gz"),
        ],
    )
    with pytest.raises(ValueError) as excinfo:
        check_shards_are_rectangular(
            conn, {"test/1/a.vcf.gz": {}, "test/1/b.vcf.gz": {}}
        )

    message = str(excinfo.value)
    assert "(2 variants)" in message
    assert "1:200:A:G" in message and "1:300:A:T" in message


def _shard(path, chromosome, samples=()):
    """A shard as check_cache_shape sees it: no zarr store behind it, and a
    chromosome that came from its position in the tree rather than its path"""
    return Shard(path=path, chromosome=chromosome, attrs={"samples": list(samples)})


def test_check_cache_shape_returns_the_chromosome(targetvariants):
    targetvariants.executemany(
        "INSERT INTO targetvariants VALUES (?, ?, ?)",
        [
            ("A", "1:100:A:C", "test/1/a.vcf.gz"),
            ("A", "1:100:A:C", "test/1/b.vcf.gz"),
        ],
    )
    shards = [
        _shard("test/1/a.vcf.gz", "1", ["s1", "s2"]),
        _shard("test/1/b.vcf.gz", "1", ["s3", "s4"]),
    ]

    chromosome = check_cache_shape(conn=targetvariants, shards=shards, sampleset="test")
    assert chromosome == "1"


def test_a_cache_holding_two_chromosomes_is_rejected(targetvariants):
    """MATCH is staged with every shard of one chromosome, so a second chromosome
    is a fan-out wiring mistake rather than a bigger job (ADR 0003)"""
    shards = [
        _shard("test/1/a.vcf.gz", "1", ["s1", "s2"]),
        _shard("test/8/a.vcf.gz", "8", ["s1", "s2"]),
    ]

    with pytest.raises(ValueError) as excinfo:
        check_cache_shape(conn=targetvariants, shards=shards, sampleset="test")

    message = str(excinfo.value)
    assert "MATCH runs per chromosome" in message
    # both chromosomes are named, so the user knows what was staged
    assert "['1', '8']" in message


def test_a_cache_with_no_shards_is_rejected(targetvariants):
    with pytest.raises(ValueError, match="No shards found under 'test'"):
        check_cache_shape(conn=targetvariants, shards=[], sampleset="test")


def _partial_scores(tmp_path, name, zarr_group, score=1.5):
    """One shard's partial_scores.parquet, as SCORE writes it"""
    path = tmp_path / name
    pl.DataFrame(
        [("test", "PGS1", 10, "HG001", 20, score, 7.0, zarr_group)],
        schema=[
            "sampleset",
            "accession",
            "n_matched",
            "sample_id",
            "allele_count",
            "score",
            "dosage_sum",
            "zarr_group",
        ],
        orient="row",
    ).write_parquet(path)
    return path


def _combine(tmp_path, paths):
    out_dir = tmp_path / "combined"
    out_dir.mkdir(exist_ok=True)
    stage = ScoreCombineStage(max_memory_gb=1, threads=1, out_dir=out_dir)
    stage.combine_scores(paths)
    with duckdb.connect(str(out_dir / "scratch" / "pgs.db"), read_only=True) as conn:
        return conn.execute(
            "SELECT sample_id, score, allele_count FROM score_table"
        ).fetchall()


def test_distinct_shards_are_summed(tmp_path):
    rows = _combine(
        tmp_path,
        [
            _partial_scores(tmp_path, "a.parquet", "test/1/a.vcf.gz", score=1.5),
            _partial_scores(tmp_path, "b.parquet", "test/8/b.vcf.gz", score=2.5),
        ],
    )
    assert rows == [("HG001", 4.0, 40)]


def test_the_same_partial_twice_is_rejected(tmp_path):
    """A glob collecting the same output twice would silently double the score"""
    path = _partial_scores(tmp_path, "a.parquet", "test/1/a.vcf.gz")

    with pytest.raises(ValueError) as excinfo:
        _combine(tmp_path, [path, path])

    message = str(excinfo.value)
    assert "test/1/a.vcf.gz (2 times)" in message
    assert "PGS1/HG001" in message


def test_two_files_holding_the_same_shard_are_rejected(tmp_path):
    """Deduplicating paths would miss this: two distinct files, one shard, which is
    what a retried task leaving a stale output looks like"""
    paths = [
        _partial_scores(tmp_path, "attempt1.parquet", "test/1/a.vcf.gz"),
        _partial_scores(tmp_path, "attempt2.parquet", "test/1/a.vcf.gz"),
    ]

    with pytest.raises(ValueError, match=r"test/1/a\.vcf\.gz \(2 times\)"):
        _combine(tmp_path, paths)
