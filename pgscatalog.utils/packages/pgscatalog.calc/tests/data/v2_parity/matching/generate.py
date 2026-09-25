"""Write matching inputs only; never replace independently obtained references."""

from __future__ import annotations

import csv
import gzip
import io
import pathlib
import tempfile

import pysam

SAMPLES = ("S0", "S1", "S2")
# Position, ID, REF, ALT, fully called diploid genotypes.
RECORDS = (
    (1000, "anchor", "A", "C", ("0/0", "0/1", "1/1")),
    (2000, "direct", "A", "C", ("0/0", "0/1", "1/1")),
    (3000, "swapped", "C", "A", ("0/1", "1/1", "0/0")),
    (4000, "flip_ref", "T", "G", ("1/1", "0/0", "0/1")),
    (5000, "flip_alt", "G", "T", ("0/0", "0/1", "1/1")),
    (6000, "no_other", "A", "C", ("0/1", "1/1", "0/0")),
    (7000, "incompatible", "C", "G", ("0/0", "0/1", "1/1")),
    (8000, "collision", "A", "C,G", ("0/1", "1/2", "2/2")),
    (9000, "distinct_alt", "A", "C,G", ("0/1", "1/2", "2/2")),
    # ALTREF is encountered first; REFALT must win, not input order.
    (11000, "candidate_alt", "C", "A", ("0/0", "0/1", "1/1")),
    (11000, "candidate_ref", "A", "C", ("1/1", "0/0", "0/1")),
    (12000, "reverse_complement", "GT", "GG", ("0/0", "0/1", "1/1")),
    (13000, "plain_complement", "TG", "GG", ("0/1", "1/1", "0/0")),
    (14000, "indel_inference", "T", "TA", ("1/1", "0/0", "0/1")),
    (15000, "forward_indel", "A", "AT", ("0/0", "0/1", "1/1")),
    (16000, "palindrome", "A", "T", ("0/0", "0/1", "1/1")),
    (17000, "no_other_ambiguous", "A", "C,G", ("0/1", "1/2", "2/2")),
)
# Position, effect allele, other allele, weight. Distinct accessions isolate policies.
SCORES = {
    "PGS990020": [
        (2000, "A", "C", 1),
        (3000, "A", "C", 2),
        (4000, "A", "C", 3),
        (5000, "A", "C", 4),
        (6000, "C", "", 5),
        (7000, "A", "C", 6),
        (999000, "A", "C", 7),
    ],
    "PGS990021": [(1000, "C", "A", 1), (8000, "A", "C", 2), (8000, "A", "G", 3)],
    "PGS990022": [(1000, "C", "A", 1), (9000, "C", "A", 2), (9000, "G", "A", 3)],
    "PGS990023": [(1000, "C", "A", 1), (11000, "A", "C", 2)],
    "PGS990024": [
        (1000, "C", "A", 1),
        (12000, "AC", "CC", 2),
        (13000, "AC", "CC", 3),
        (14000, "A", "AT", 4),
        (15000, "AT", "A", 5),
    ],
    "PGS990025": [(1000, "C", "A", 1), (16000, "A", "T", 2), (17000, "C", "", 3)],
}


def generate(out_dir: pathlib.Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    header = (
        "##fileformat=VCFv4.2\n##contig=<ID=1>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(SAMPLES)
        + "\n"
    )
    body = "".join(
        f"1\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gt) + "\n"
        for pos, vid, ref, alt, gt in RECORDS
    )
    with tempfile.TemporaryDirectory() as temporary:
        plain = pathlib.Path(temporary) / "target.vcf"
        plain.write_text(header + body)
        target = out_dir / "target.chr1.vcf.gz"
        pysam.tabix_compress(str(plain), str(target), force=True)
        pysam.tabix_index(str(target), preset="vcf", force=True)
    for accession, variants in SCORES.items():
        columns = [
            "chr_name",
            "chr_position",
            "effect_allele",
            "other_allele",
            "effect_weight",
        ]
        custom = io.StringIO()
        writer = csv.writer(custom, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(("1", *variant) for variant in variants)
        (out_dir / f"{accession}.txt").write_text(
            f"#pgs_id={accession}\n#pgs_name=MatchingParity\n"
            "#trait_reported=Synthetic test fixture\n#genome_build=GRCh38\n"
            + custom.getvalue()
        )
        processed = io.StringIO()
        writer = csv.writer(processed, delimiter="\t", lineterminator="\n")
        writer.writerow(
            columns + ["effect_type", "is_duplicated", "accession", "row_nr"]
        )
        writer.writerows(
            ("1", *variant, "additive", "False", accession, row)
            for row, variant in enumerate(variants)
        )
        with (
            (out_dir / f"{accession}.txt.gz").open("wb") as raw,
            gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle,
        ):
            handle.write(processed.getvalue().encode())


if __name__ == "__main__":
    generate(pathlib.Path(__file__).parent)
