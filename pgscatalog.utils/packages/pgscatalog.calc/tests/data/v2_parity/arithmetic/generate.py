"""Write deterministic hard-call arithmetic inputs; leave references untouched."""

import argparse
import gzip
import tempfile
from pathlib import Path

import pysam

SAMPLES = ("S0", "S1", "S2", "S3")
# Equal calls at the first pair make large opposite weights cancel exactly.
CALLS = (
    (0, 1, 2, 1),
    (0, 1, 2, 1),
    (2, 0, 1, 2),
    (1, 2, 0, 0),
    (2, 1, 0, 1),
    (0, 2, 1, 0),
)
# position index, effect allele, decimal weight, genetic effect model.
ROWS = {
    "PGS990001": (
        (1, "C", "123.4567", "additive"),
        (2, "C", "-123.4567", "additive"),
        (3, "A", "0.000123", "additive"),
        (4, "C", "-0.17", "additive"),
    ),
    "PGS990002": (
        (1, "A", "-0.31", "additive"),
        (3, "A", "0.27", "additive"),
        (5, "C", "1.234567", "additive"),
    ),
    "PGS990003": ((6, "C", "-0.000317", "additive"),),
    "PGS990004": (
        (1, "C", "0.13", "additive"),
        (2, "A", "-0.29", "dominant"),
        (3, "C", "0.37", "recessive"),
    ),
}


def generate(out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    for chrom in ("1", "2"):
        lines = [
            "##fileformat=VCFv4.2\n",
            f"##contig=<ID={chrom}>\n",
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(SAMPLES)
            + "\n",
        ]
        for index, counts in enumerate(CALLS, 1):
            # Different chromosome contributions expose compensating shard errors.
            counts = counts if chrom == "1" else counts[::-1]
            calls = [("0/0", "0/1", "1/1")[x] for x in counts]
            lines.append(
                f"{chrom}\t{index * 1000}\tarith_{chrom}_{index}"
                "\tA\tC\t.\tPASS\t.\tGT\t" + "\t".join(calls) + "\n"
            )
        target = out_dir / f"target.chr{chrom}.vcf.gz"
        with tempfile.TemporaryDirectory() as temp:
            plain = Path(temp) / "target.vcf"
            plain.write_text("".join(lines))
            pysam.tabix_compress(str(plain), str(target), force=True)
        pysam.tabix_index(str(target), preset="vcf", force=True)
    for accession, rows in ROWS.items():
        custom = [
            f"#pgs_id={accession}\n#pgs_name=ArithmeticParity\n"
            "#trait_reported=Synthetic test fixture\n#genome_build=GRCh38\n",
            "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight"
            "\tis_dominant\tis_recessive\n",
        ]
        processed = [
            "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight"
            "\teffect_type\tis_duplicated\taccession\trow_nr\n"
        ]
        for chrom in ("1", "2"):
            for index, effect, weight, model in rows:
                other = "A" if effect == "C" else "C"
                fields = f"{chrom}\t{index * 1000}\t{effect}\t{other}\t{weight}"
                custom.append(
                    fields + f"\t{model == 'dominant'}\t{model == 'recessive'}\n"
                )
                processed.append(
                    fields + f"\t{model}\tFalse\t{accession}\t{len(processed) - 1}\n"
                )
        (out_dir / f"{accession}.txt").write_text("".join(custom))
        with (
            (out_dir / f"{accession}.txt.gz").open("wb") as raw,
            gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle,
        ):
            handle.write("".join(processed).encode())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out_dir", type=Path)
    generate(parser.parse_args().out_dir)
