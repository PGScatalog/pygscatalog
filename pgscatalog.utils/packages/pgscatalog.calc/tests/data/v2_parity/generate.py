"""Generate synthetic scoring inputs; never generate the PLINK reference result."""

from __future__ import annotations

import argparse
import gzip
import pathlib
import tempfile

import pysam

CHROMS = (*map(str, range(1, 23)), "X", "Y", "MT")
SAMPLES = ("XX_1", "XX_2", "XY_1", "XY_2")
ACCESSION = "PGS999997"
# kind, REF, ALT, diploid calls, haploid calls (in SAMPLES order)
VARIANTS = (
    ("snp", "A", "C", ("0/0", "0/1", "1/1", "0/1"), ("0", "1", "1", "0")),
    ("multi", "A", "C,G", ("0/0", "0/1", "1/2", "2/2"), ("0", "1", "1", "2")),
    ("ins", "A", "AT", ("0/1", "1/1", "0/0", "0/1"), ()),
    ("del", "AC", "A", ("1/1", "0/0", "0/1", "0/0"), ()),
)


def generate(out_dir: pathlib.Path) -> None:
    """Write fixed inputs, reproducible within the same dependency environment."""
    out_dir.mkdir(parents=True, exist_ok=True)
    scores = [
        "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\t"
        "effect_type\tis_duplicated\taccession\trow_nr\n"
    ]
    weights = ["ID\tALLELE\tWEIGHT\n"]
    for chrom_index, chrom in enumerate(CHROMS):
        samples = SAMPLES[2:] if chrom == "Y" else SAMPLES
        lines = [
            "##fileformat=VCFv4.2\n",
            "##source=pygscatalog_synthetic_parity\n",
            f"##contig=<ID={chrom}>\n",
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(samples)
            + "\n",
        ]
        for index, (kind, ref, alt, diploid, haploid) in enumerate(VARIANTS):
            if index >= 2 and chrom != "1":
                continue
            # X/Y are well outside GRCh38 PARs; MT positions fit its length.
            pos = (10_000_000 if chrom in ("X", "Y") else 0) + (index + 1) * 1000
            variant_id = f"syn_{chrom}_{kind}"
            alleles = [ref, *alt.split(",")]
            effect_index = chrom_index % len(alleles) if index < 2 else 1
            effect = alleles[effect_index]
            other = alleles[1] if effect_index == 0 else ref
            # Binary-exact fractions keep this fixture focused on allele counting.
            weight = (chrom_index + 1) * (index + 1) / 16
            if index % 2:
                weight = -weight
            calls = diploid
            if chrom == "X":
                calls = (*diploid[:2], *haploid[2:])
            elif chrom == "Y":
                calls = haploid[2:]
            elif chrom == "MT":
                calls = haploid
            lines.append(
                f"{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t"
                + "\t".join(calls)
                + "\n"
            )
            scores.append(
                f"{chrom}\t{pos}\t{effect}\t{other}\t{weight}\tadditive\tFalse\t"
                f"{ACCESSION}\t{len(scores) - 1}\n"
            )
            weights.append(f"{variant_id}\t{effect}\t{weight}\n")
        target = out_dir / f"target.chr{chrom}.vcf.gz"
        with tempfile.TemporaryDirectory() as temp:
            plain = pathlib.Path(temp) / "target.vcf"
            plain.write_text("".join(lines), encoding="ascii")
            pysam.tabix_compress(str(plain), str(target), force=True)
        pysam.tabix_index(str(target), preset="vcf", force=True)

    with (
        (out_dir / "score.txt.gz").open("wb") as raw,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed,
    ):
        compressed.write("".join(scores).encode("ascii"))
    (out_dir / "plink.weights.tsv").write_text("".join(weights), encoding="ascii")
    declarations = [(sample, 2 if sample.startswith("XX") else 1) for sample in SAMPLES]
    rows = [f"{sample}\t{sex}\n" for sample, sex in declarations]
    (out_dir / "sex.tsv").write_text(
        "sample_id\tsex\n" + "".join(rows), encoding="ascii"
    )
    (out_dir / "samples.psam").write_text(
        "#IID\tSEX\n" + "".join(rows), encoding="ascii"
    )
    (out_dir / "samples.Y.psam").write_text(
        "#IID\tSEX\n" + "".join(rows[2:]), encoding="ascii"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "out_dir", type=pathlib.Path, help="Directory for fixture inputs"
    )
    generate(parser.parse_args().out_dir)
