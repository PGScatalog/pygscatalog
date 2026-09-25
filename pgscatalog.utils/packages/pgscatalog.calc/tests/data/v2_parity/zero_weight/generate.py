"""Generate zero-weight score inputs; external references are never overwritten."""

import argparse
import gzip
import tempfile
from pathlib import Path

import pysam


def generate(out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    content = (
        "##fileformat=VCFv4.2\n##contig=<ID=1>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\tS1\tS2\n"
        "1\t1000\tzero_1\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n"
        "1\t2000\tzero_2\tA\tC\t.\tPASS\t.\tGT\t1/1\t0/1\t0/0\n"
    )
    target = out_dir / "target.chr1.vcf.gz"
    with tempfile.TemporaryDirectory() as temp:
        plain = Path(temp) / "target.vcf"
        plain.write_text(content, encoding="ascii")
        pysam.tabix_compress(str(plain), str(target), force=True)
    pysam.tabix_index(str(target), preset="vcf", force=True)
    header = "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight"
    for accession, weights in (("PGS990040", (0, 2)), ("PGS990041", (0, 0))):
        custom = (
            f"#pgs_id={accession}\n#pgs_name=ZeroWeights\n"
            "#trait_reported=Synthetic test fixture\n#genome_build=GRCh38\n"
            + header
            + "\n"
        )
        processed = header + "\teffect_type\tis_duplicated\taccession\trow_nr\n"
        for index, weight in enumerate(weights):
            row = f"1\t{(index + 1) * 1000}\tC\tA\t{weight}"
            custom += row + "\n"
            processed += row + f"\tadditive\tFalse\t{accession}\t{index}\n"
        (out_dir / f"{accession}.txt").write_text(custom, encoding="ascii")
        with (
            (out_dir / f"{accession}.txt.gz").open("wb") as raw,
            gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as stream,
        ):
            stream.write(processed.encode("ascii"))
    (out_dir / "samples.psam").write_text(
        "#IID\tSEX\nS0\t2\nS1\t2\nS2\t2\n", encoding="ascii"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out_dir", type=Path)
    generate(parser.parse_args().out_dir)
