"""Acceptance tests for haploid calls, end to end: chrX from VCF (ADR 0007).

The fixture VCF already holds real 1000G chrX data - 1599 men whose
non-pseudoautosomal calls are written `1`, and 1603 women whose calls are written
`0/1` - so ploidy is observed from the file at every layer, from the ragged `GT`
to a per-sample allele count in the output. Nothing here knows where the
pseudoautosomal region is.
"""

import gzip
import pathlib

import numpy as np
import polars as pl
import pysam
import pytest
import zarr
import zarr.storage

from pgscatalog.calc.cli.load import unzip_zarr
from pgscatalog.calc.lib.constants import MISSING_GENOTYPE_SENTINEL_VALUE

from .conftest import SAMPLESET, _is_ambiguous, _write_vcf, load_target, run_serial

# GRCh38 PAR1 ends here, which is only used to place a variant the caller wrote as
# diploid for everyone. Nothing in the pipeline reads it
PAR_POSITION = 2781479

# a position the target doesn't hold, so LOAD fabricates a padding row for it
ABSENT_POSITION = 99_999_999


def _write_scorefile(path: pathlib.Path, records: list, extra: list = ()) -> None:
    """One scoring file variant per record, effect allele = ALT, weight = 1.

    A weight of 1 makes score the sum of effect allele dosages, so the numbers the
    acceptance criteria talk about are readable in the output.
    """
    with gzip.open(path, "wt") as f:
        f.write(
            "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\t"
            "effect_type\tis_duplicated\taccession\trow_nr\n"
        )
        rows = [(rec.pos, rec.alts[0], rec.ref) for rec in records] + list(extra)
        for row_nr, (pos, effect, other) in enumerate(rows):
            f.write(
                f"X\t{pos}\t{effect}\t{other}\t1.0\tadditive\tFalse\tPGS999995\t{row_nr}\n"
            )


def _observed(records: list) -> tuple[dict[str, int], dict[str, int]]:
    """What the file says: each sample's ALT dosage and allele count"""
    dosage: dict[str, int] = {}
    alleles: dict[str, int] = {}
    for rec in records:
        for name, call in rec.samples.items():
            gt = call["GT"]
            dosage[name] = dosage.get(name, 0) + sum(allele == 1 for allele in gt)
            alleles[name] = alleles.get(name, 0) + len(gt)
    return dosage, alleles


@pytest.fixture(scope="session")
def chrx_records(vcf_path):
    """Simple biallelic chrX SNPs, plus one the caller wrote as diploid for every
    sample (a pseudoautosomal variant)"""
    with pysam.VariantFile(str(vcf_path)) as vcf:
        header = vcf.header.copy()
        records = [
            rec
            for rec in vcf.fetch("X")
            if len(rec.ref) == 1
            and len(rec.alts) == 1
            and len(rec.alts[0]) == 1
            and not _is_ambiguous(rec.ref, rec.alts[0])
        ]

    par, records = records[0], records[1:]
    par.pos = PAR_POSITION
    for call in par.samples.values():
        if len(call["GT"]) == 1:
            call["GT"] = (call["GT"][0], 0)

    assert len(records) > 10, "not enough chrX test variants"
    return header, par, records


@pytest.fixture(scope="session")
def chrx_target(chrx_records, tmp_path_factory) -> pathlib.Path:
    header, par, records = chrx_records
    path = tmp_path_factory.mktemp("chrx") / "chrX.vcf.gz"
    _write_vcf(path, header, [par, *records])
    return path


@pytest.fixture(scope="session")
def chrx_scorefile(chrx_records, tmp_path_factory) -> pathlib.Path:
    """Covers every variant in the target, plus one position it doesn't hold"""
    _, par, records = chrx_records
    path = tmp_path_factory.mktemp("chrx_score") / "PGS999995.txt.gz"
    _write_scorefile(path, [par, *records], extra=[(ABSENT_POSITION, "A", "G")])
    return path


def test_mixed_sex_chrx_scores_end_to_end(
    chrx_records, chrx_target, chrx_scorefile, sex_map, tmp_path, monkeypatch
):
    """A hemizygous effect allele contributes 1, a female heterozygote 1 and a
    female homozygote 2, and each sample's allele count follows its own ploidy.

    The padded position is in here too: LOAD fabricates a diploid-shaped
    all-missing row for it, and if anything read ploidy from that row every
    sample's allele count would be two too high.
    """
    _, par, records = chrx_records
    expected_dosage, expected_alleles = _observed([par, *records])

    scores = run_serial(
        [(chrx_target, "X")], tmp_path, chrx_scorefile, monkeypatch, sex_map=sex_map
    )

    assert scores.shape[0] == 3202
    # the padded position can't match, so it counts against the match rate and
    # contributes no alleles to anybody
    assert scores["n_matched"].unique().to_list() == [len(records) + 1]

    for row in scores.iter_rows(named=True):
        sample = row["sample_id"]
        assert row["dosage_sum"] == expected_dosage[sample], sample
        assert row["allele_count"] == expected_alleles[sample], sample

    # a weight of 1 per variant, so the score is the dosage sum and the average
    # weighted score divides it by the sample's own allele count
    assert scores["score"].to_list() == pytest.approx(scores["dosage_sum"].to_list())
    # scores are written with 6 decimal places to limit identifiability
    assert scores["score_avg"].to_list() == pytest.approx(
        (scores["score"] / scores["allele_count"]).to_list(), abs=1e-6
    )

    # men carry one allele per non-pseudoautosomal variant and two at the
    # pseudoautosomal one; women carry two everywhere
    hemizygous = {s for s, n in expected_alleles.items() if n < 2 * (len(records) + 1)}
    assert len(hemizygous) == 1599
    assert {expected_alleles[s] for s in hemizygous} == {len(records) + 2}
    assert scores.filter(pl.col("sample_id").is_in(hemizygous))["dosage_sum"].max() <= (
        len(records) + 2
    )


def test_a_pseudoautosomal_variant_is_diploid_for_every_sample(
    chrx_records, chrx_target, sex_map, tmp_path, monkeypatch
):
    """No PAR coordinate map is consulted: the same men are hemizygous a few
    kilobases away, and diploid here, because that is what the file says"""
    _, par, _ = chrx_records
    scorefile = tmp_path / "PGS999994.txt.gz"
    _write_scorefile(scorefile, [par])

    scores = run_serial(
        [(chrx_target, "X")], tmp_path, scorefile, monkeypatch, sex_map=sex_map
    )

    assert scores["allele_count"].unique().to_list() == [2]
    assert scores["dosage_sum"].max() == 2


def test_padded_positions_are_diploid_shaped_but_carry_no_alleles(
    chrx_target, chrx_scorefile, sex_map, tmp_path, monkeypatch
):
    """The other half of the padding evidence: the fabricated row really is width
    2 on a chromosome where most calls are hemizygous. It's harmless because it
    holds no alleles and has no REF, so it never becomes a match candidate"""
    zip_path = load_target(
        f"{chrx_target}:X",
        tmp_path / "cache",
        chrx_scorefile,
        monkeypatch,
        sex_map=sex_map,
    )
    unzip_zarr(zip_path, tmp_path / "unzipped")
    leaf = f"{SAMPLESET}/X/{chrx_target.name}"
    store = zarr.storage.LocalStore(tmp_path / "unzipped" / "genotypes.zarr")
    root = zarr.open_group(store=store, mode="r")

    positions = root[f"{leaf}/meta/chr_pos"][:]
    refs = root[f"{leaf}/meta/ref"][:]
    padded = int(np.flatnonzero(positions == ABSENT_POSITION)[0])
    genotypes = root[f"{leaf}/genotypes"][padded]

    assert refs[padded] == ""  # no REF, so update_match_table skips the position
    assert genotypes.shape[1] == 2
    assert np.all(genotypes == MISSING_GENOTYPE_SENTINEL_VALUE)
