"""The fixture stays reproducible and agrees with its two scoring inputs."""

from __future__ import annotations

import csv
import gzip
import pathlib
import runpy

import pysam


def test_generate(tmp_path):
    generate = runpy.run_path(str(pathlib.Path(__file__).with_name("generate.py")))[
        "generate"
    ]
    first, second = tmp_path / "first", tmp_path / "second"
    generate(first)
    generate(second)
    assert len(list(first.iterdir())) == 55
    assert {p.name: p.read_bytes() for p in first.iterdir()} == {
        p.name: p.read_bytes() for p in second.iterdir()
    }
    # Regeneration must leave the independently obtained reference untouched.
    expected = first / "expected.csv"
    expected.write_text("independent reference\n")
    generate(first)
    assert expected.read_text() == "independent reference\n"

    with gzip.open(first / "score.txt.gz", "rt") as handle:
        scores = list(csv.DictReader(handle, delimiter="\t"))
    with (first / "plink.weights.tsv").open() as handle:
        weights = list(csv.DictReader(handle, delimiter="\t"))
    assert len(scores) == len(weights) == 52
    by_position = {(r["chr_name"], int(r["chr_position"])): r for r in scores}
    by_id = {r["ID"]: r for r in weights}
    assert len(by_position) == len(by_id) == 52
    assert [int(r["row_nr"]) for r in scores] == list(range(52))
    assert {r["accession"] for r in scores} == {"PGS999997"}
    assert {r["effect_type"] for r in scores} == {"additive"}
    observed_ids, effect_indices, dosages = set(), set(), set()
    chroms = [*map(str, range(1, 23)), "X", "Y", "MT"]
    samples = ["XX_1", "XX_2", "XY_1", "XY_2"]
    assert (first / "sex.tsv").read_text() == (
        "sample_id\tsex\nXX_1\t2\nXX_2\t2\nXY_1\t1\nXY_2\t1\n"
    )
    assert (first / "samples.psam").read_text() == (
        "#IID\tSEX\nXX_1\t2\nXX_2\t2\nXY_1\t1\nXY_2\t1\n"
    )
    assert (first / "samples.Y.psam").read_text() == "#IID\tSEX\nXY_1\t1\nXY_2\t1\n"
    for chrom in chroms:
        with pysam.VariantFile(first / f"target.chr{chrom}.vcf.gz") as vcf:
            assert list(vcf.header.samples) == (
                samples[2:] if chrom == "Y" else samples
            )
            records = list(vcf.fetch(chrom))  # Requires a working index.
            assert len(records) == (4 if chrom == "1" else 2)
            assert [r.pos for r in records] == sorted(r.pos for r in records)
            assert [len(r.alleles) for r in records[:2]] == [2, 3]
            for record in records:
                observed_ids.add(record.id)
                score = by_position[chrom, record.pos]
                weight = by_id[record.id]
                assert score["effect_allele"] == weight["ALLELE"]
                assert score["effect_weight"] == weight["WEIGHT"]
                assert score["other_allele"] in record.alleles
                effect_index = record.alleles.index(weight["ALLELE"])
                effect_indices.add(effect_index)
                assert all(set(a) <= set("ACTG") for a in record.alleles)
                for sample, call in record.samples.items():
                    gt = call["GT"]
                    haploid = chrom in ("Y", "MT") or (
                        chrom == "X" and sample.startswith("XY")
                    )
                    assert len(gt) == (1 if haploid else 2)
                    assert all(
                        a is not None and 0 <= a < len(record.alleles) for a in gt
                    )
                    dosages.add(gt.count(effect_index))
            if chrom == "1":
                assert len(records[2].ref) < len(records[2].alts[0])
                assert len(records[3].ref) > len(records[3].alts[0])
                assert records[1].samples["XY_1"]["GT"] == (1, 2)
    assert observed_ids == set(by_id)
    assert effect_indices == dosages == {0, 1, 2}
    assert any(float(r["WEIGHT"]) < 0 for r in weights)
    assert any(float(r["WEIGHT"]) > 0 for r in weights)
