"""Technical sample partitions must not change imputation (ADR 0016)."""

import gzip

import numpy as np
import polars as pl
import pysam
import pytest

from .conftest import _subset_samples, _write_vcf, run_serial


def _inputs(tmp_path, n_samples=60, all_missing=False):
    header = pysam.VariantHeader()
    for chrom in ("1", "2"):
        header.contigs.add(chrom)
    header.formats.add("GT", 1, "String", "Genotype")
    for i in range(n_samples):
        header.add_sample(f"sample{i}")
    records = []
    for chrom, alt in (("1", "G"), ("2", "C")):
        record = header.new_record(contig=chrom, start=9, alleles=("A", alt))
        for i, sample in enumerate(header.samples):
            if all_missing or i < 10:
                gt = (None, None) if i else (0, None)
            else:
                allele = int(i >= 20) if chrom == "1" else int(i < 20)
                gt = (allele, allele)
            record.samples[sample]["GT"] = gt
        records.append(record)
    target = tmp_path / "whole.vcf.gz"
    _write_vcf(target, header, records)
    scorefile = tmp_path / "PGS999998.txt.gz"
    with gzip.open(scorefile, "wt") as out:
        out.write(
            "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\t"
            "effect_type\tis_duplicated\taccession\trow_nr\n"
        )
        for i, (chrom, alt) in enumerate((("1", "G"), ("2", "C"))):
            out.write(
                f"{chrom}\t10\t{alt}\tA\t{i + 1}\tadditive\tFalse\tPGS999998\t{i}\n"
            )
    return target, scorefile, header, records


@pytest.mark.parametrize("by_chromosome", [False, True])
def test_sample_partitions_preserve_scores_and_dosages(
    tmp_path, monkeypatch, by_chromosome
):
    target, scorefile, header, records = _inputs(tmp_path)
    whole = run_serial([target], tmp_path / "whole", scorefile, monkeypatch)
    samples = list(header.samples)
    targets = []
    sources = [(target, None)]
    if by_chromosome:
        sources = []
        for record in records:
            path = tmp_path / f"chr{record.chrom}.vcf.gz"
            _write_vcf(path, header, [record])
            sources.append((path, record.chrom))
    for source, chrom in sources:
        # Unequal batches; the first has no complete calls at either variant.
        for i, batch in enumerate((samples[:10], samples[10:20], samples[20:])):
            path = tmp_path / f"{chrom}.batch{i}.vcf.gz"
            _subset_samples(source, path, batch)
            targets.append((path, chrom) if chrom else path)
    split = run_serial(targets, tmp_path / "split", scorefile, monkeypatch)
    whole, split = whole.sort("sample_id"), split.sort("sample_id")
    assert whole["sample_id"].to_list() == split["sample_id"].to_list()
    for column in ("score", "dosage_sum", "allele_count", "n_matched"):
        np.testing.assert_allclose(whole[column], split[column], rtol=0, atol=1e-12)
    # Frequencies are 40/50 and 10/50: weights 1 and 2 give 1.6 + 0.8.
    assert split.filter(pl.col("sample_id") == "sample0")["score"][0] == pytest.approx(
        2.4
    )
    matches = pl.read_parquet(
        tmp_path / "split/stages/match_combine/match_table.parquet"
    )
    assert matches["sampleset_n_samples"].unique().to_list() == [60]
    for row_nr, freq in ((0, 0.8), (1, 0.2)):
        assert (
            matches.filter(pl.col("row_nr") == row_nr)[
                "effect_allele_frequency"
            ].to_list()
            == [freq] * 3
        )


def test_floor_counts_distinct_samples_across_chromosomes(tmp_path, monkeypatch):
    _, scorefile, header, records = _inputs(tmp_path, n_samples=49)
    targets = []
    for record in records:
        path = tmp_path / f"chr{record.chrom}.vcf.gz"
        _write_vcf(path, header, [record])
        targets.append((path, record.chrom))
    with pytest.raises(
        ValueError, match="n_samples=49 is less than n_minimum_samples=50"
    ):
        run_serial(targets, tmp_path / "run", scorefile, monkeypatch)


def test_no_called_alleles_across_sampleset_names_variant(tmp_path, monkeypatch):
    target, scorefile, _, _ = _inputs(tmp_path, all_missing=True)
    with pytest.raises(
        ValueError, match="no called alleles across the sampleset.*1:10:A:G"
    ):
        run_serial([target], tmp_path / "run", scorefile, monkeypatch)
