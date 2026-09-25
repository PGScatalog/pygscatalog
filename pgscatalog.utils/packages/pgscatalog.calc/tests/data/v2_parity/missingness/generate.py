"""Generate hard-call missingness inputs, never external reference outputs."""

import argparse
import gzip
import tempfile
from pathlib import Path

import pysam


def generate(out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    samples = [f"S{i:02d}" for i in range(52)]
    for name, chrom, selected, missing in (
        ("target.chr1", "1", samples, {"S00", "S01"}),
        ("target.chr2", "2", samples, {"S00", "S02"}),
        ("all_uncalled.chr1", "1", samples, set(samples)),
        ("too_small.chr1", "1", samples[:4], {"S00", "S01"}),
    ):
        content = (
            "##fileformat=VCFv4.2\n"
            f"##contig=<ID={chrom}>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(selected)
            + "\n"
        )
        for pos, called in ((1000, "0/1"), (2000, "0/0")):
            content += (
                f"{chrom}\t{pos}\tmissing_{chrom}_{pos}\tA\tC\t.\tPASS\t.\tGT\t"
                + "\t".join(
                    "./." if sample in missing else called for sample in selected
                )
                + "\n"
            )
        target = out_dir / f"{name}.vcf.gz"
        with tempfile.TemporaryDirectory() as temp:
            plain = Path(temp) / "target.vcf"
            plain.write_text(content, encoding="ascii")
            pysam.tabix_compress(str(plain), str(target), force=True)
        pysam.tabix_index(str(target), preset="vcf", force=True)

    rows = [
        ("1", 1000, "C", "A", "2"),
        ("1", 2000, "C", "A", "0.25"),
        ("2", 1000, "A", "C", "-0.5"),
        ("2", 2000, "A", "C", "0.5"),
    ]
    custom = (
        "#pgs_id=PGS990010\n#pgs_name=IntegerMissingness\n"
        "#trait_reported=Synthetic test fixture\n#genome_build=GRCh38\n"
        "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n"
    )
    processed = (
        "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\t"
        "effect_type\tis_duplicated\taccession\trow_nr\n"
    )
    for i, (chrom, position, effect, other, weight) in enumerate(rows):
        row = f"{chrom}\t{position}\t{effect}\t{other}\t{weight}"
        custom += row + "\n"
        processed += row + f"\tadditive\tFalse\tPGS990010\t{i}\n"
    (out_dir / "PGS990010.txt").write_text(custom, encoding="ascii")
    with (
        (out_dir / "score.txt.gz").open("wb") as raw,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as stream,
    ):
        stream.write(processed.encode("ascii"))
    (out_dir / "samples.psam").write_text(
        "#IID\tSEX\n" + "".join(f"{sample}\t2\n" for sample in samples),
        encoding="ascii",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out_dir", type=Path)
    generate(parser.parse_args().out_dir)
