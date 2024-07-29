import argparse
import logging
import pathlib

from xopen import xopen
import csv
import textwrap
import heapq
import tempfile


logger = logging.getLogger(__name__)


def run_intersect():
    args = parse_args()
    outdir = pathlib.Path(args.outdir)

    if not outdir.is_dir():
        raise NotADirectoryError

    tmpdir = tempfile.mkdtemp()

    if args.verbose:
        logging.getLogger("pgscatalog.core").setLevel(logging.DEBUG)
        logging.getLogger("pgscatalog.match").setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    # Process & sort reference variants
    count_var_r = 0
    o_tmp_r = []
    logger.info("Reading REFERENCE variants: {}".format(args.reference))
    ref_heap = []
    ref_pvar = read_var_general(args.reference, chrom=args.filter_chrom)
    for v in ref_pvar:
        count_var_r += 1
        ALTs = v["ALT"].split(",")
        IS_MA_REF = len(ALTs) > 1
        for i, ALT in enumerate(ALTs):
            if v["REF"] < ALT:
                key = "{}:{}:{}:{}".format(v["#CHROM"], v["POS"], v["REF"], ALT)
            else:
                key = "{}:{}:{}:{}".format(v["#CHROM"], v["POS"], ALT, v["REF"])

            IS_INDEL = (len(v["REF"]) > 1) | (len(ALT) > 1)
            STRANDAMB = v["REF"] == allele_complement(ALT)
            ref_heap.append(
                ([key, v["ID"], v["REF"]], [IS_INDEL, STRANDAMB, IS_MA_REF])
            )

        if count_var_r % 500000 == 0:
            heapq.heapify(ref_heap)
            tmppath = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)
            with open(tmppath.name, "wt") as outf:
                o_tmp_r.append(tmppath.name)
                outf.write(
                    "CHR:POS:A0:A1\tID_REF\tREF_REF\tIS_INDEL\tSTRANDAMB\tIS_MA_REF\n"
                )
                for i in range(len(ref_heap)):
                    popped = heapq.heappop(ref_heap)
                    outf.write(
                        "\t".join([str(x) for x in popped[0] + popped[1]]) + "\n"
                    )
            ref_heap = []
            logger.info("Processed {} REFERENCE variants".format(count_var_r))

    if len(ref_heap) > 0:
        heapq.heapify(ref_heap)
        tmppath = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)
        with open(tmppath.name, "wt") as outf:
            o_tmp_r.append(tmppath.name)
            outf.write(
                "CHR:POS:A0:A1\tID_REF\tREF_REF\tIS_INDEL\tSTRANDAMB\tIS_MA_REF\n"
            )
            for i in range(len(ref_heap)):
                popped = heapq.heappop(ref_heap)
                outf.write("\t".join([str(x) for x in popped[0] + popped[1]]) + "\n")
        del ref_heap
        logger.info("Processed {} REFERENCE variants".format(count_var_r))

    logger.info("Outputting REFERNCE variants -> reference_variants.txt.gz")
    with xopen(outdir / "reference_variants.txt.gz", "wt") as outf:
        outf.write("CHR:POS:A0:A1\tID_REF\tREF_REF\tIS_INDEL\tSTRANDAMB\tIS_MA_REF\n")
        for v in heapq.merge(
            *[read_var_general(x) for x in o_tmp_r],
            key=lambda v: (v["CHR:POS:A0:A1"], v["ID_REF"], v["REF_REF"]),
        ):
            outf.write("\t".join(v.values()) + "\n")

    # Process & sort target variants
    count_var_t = 0
    o_tmp_t = []
    target_heap = []
    for path in args.target:
        # assume inputs always have two extensions (bim.zst / pvar.zst)
        path = pathlib.Path(path)
        basename = path.name
        stem = pathlib.Path(basename.split(".")[0])

        logger.info("Reading TARGET variants: {}".format(path))
        pvar = read_var_general(
            path, chrom=None
        )  # essential not to filter target (messes up common line indexing)

        loc_afreq = path.parent / stem.with_suffix(".afreq.gz")
        afreq = read_var_general(
            loc_afreq, chrom=None
        )  # essential not to filter target (messes up common line indexing)

        loc_vmiss = path.parent / stem.with_suffix(".vmiss.gz")
        vmiss = read_var_general(
            loc_vmiss, chrom=None
        )  # essential not to filter target (messes up common line indexing)

        for v, freq, miss in zip(pvar, afreq, vmiss):
            if not all([v["ID"], freq["ID"], miss["#ID"]]):
                raise ValueError("TARGET variant files are not sorted")

            count_var_t += 1
            ALTs = v["ALT"].split(",")
            ALT_FREQS = [float(x) for x in freq["ALT_FREQS"].split(",")]
            F_MISS_DOSAGE = miss["F_MISS_DOSAGE"]
            IS_MA_TARGET = len(ALTs) > 1
            for i, ALT in enumerate(ALTs):
                if v["REF"] < ALT:
                    key = "{}:{}:{}:{}".format(v["#CHROM"], v["POS"], v["REF"], ALT)
                else:
                    key = "{}:{}:{}:{}".format(v["#CHROM"], v["POS"], ALT, v["REF"])
                target_heap.append(
                    (
                        [key, v["ID"], v["REF"]],
                        [IS_MA_TARGET, ALT_FREQS[i], F_MISS_DOSAGE],
                    )
                )

            if count_var_t % 500000 == 0:
                heapq.heapify(target_heap)
                tmppath = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)
                with open(tmppath.name, "wt") as outf:
                    o_tmp_t.append(tmppath.name)
                    outf.write(
                        "CHR:POS:A0:A1\tID_TARGET\tREF_TARGET\tIS_MA_TARGET\tAAF\tF_MISS_DOSAGE\n"
                    )
                    for i in range(len(target_heap)):
                        popped = heapq.heappop(target_heap)
                        outf.write(
                            "\t".join([str(x) for x in popped[0] + popped[1]]) + "\n"
                        )
                target_heap = []
                logger.info("Processed {} TARGET variants".format(count_var_t))

    if len(target_heap) > 0:
        heapq.heapify(target_heap)
        tmppath = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)
        with open(tmppath.name, "wt") as outf:
            o_tmp_t.append(tmppath.name)
            outf.write(
                "CHR:POS:A0:A1\tID_TARGET\tREF_TARGET\tIS_MA_TARGET\tAAF\tF_MISS_DOSAGE\n"
            )
            for i in range(len(target_heap)):
                popped = heapq.heappop(target_heap)
                outf.write("\t".join([str(x) for x in popped[0] + popped[1]]) + "\n")
        del target_heap
        logger.info("Processed {} TARGET variants".format(count_var_t))

    logger.info("Outputting TARGET variants -> target_variants.txt.gz")
    with xopen(outdir / "target_variants.txt.gz", "wt") as outf:
        outf.write(
            "CHR:POS:A0:A1\tID_TARGET\tREF_TARGET\tIS_MA_TARGET\tAAF\tF_MISS_DOSAGE\n"
        )
        for v in heapq.merge(
            *[read_var_general(x) for x in o_tmp_t],
            key=lambda v: (v["CHR:POS:A0:A1"], v["ID_TARGET"], v["REF_TARGET"]),
        ):
            outf.write("\t".join(v.values()) + "\n")

    # Merge matched variants on sorted files
    logger.info("Joining & outputting matched variants -> matched_variants.txt.gz")
    n_matched = 0
    n_PCA_ELIGIBLE = 0
    with xopen(outdir / "matched_variants.txt.gz", "w") as csvfile:
        for vmatch in sorted_join_variants(
            outdir / "reference_variants.txt.gz", outdir / "target_variants.txt.gz"
        ):
            n_matched += 1
            vmatch["SAME_REF"] = vmatch["REF_REF"] == vmatch["REF_REF"]

            # Define variant's eligibility for PCA
            # From original implementation: ((IS_MA_REF == FALSE) && (IS_MA_TARGET == FALSE)) && (((IS_INDEL == FALSE) && (STRANDAMB == FALSE)) || ((IS_INDEL == TRUE) && (SAME_REF == TRUE)))
            PCA_ELIGIBLE = (
                (vmatch["IS_MA_REF"] == "False") and (vmatch["IS_MA_TARGET"] == "False")
            ) and (
                ((vmatch["IS_INDEL"] == "False") and (vmatch["STRANDAMB"] == "False"))
                or ((vmatch["IS_INDEL"] == "True") and (vmatch["SAME_REF"] == "True"))
            )

            PCA_ELIGIBLE = (
                PCA_ELIGIBLE
                and (aaf2maf(float(vmatch["AAF"])) >= args.maf_filter)
                and (float(vmatch["F_MISS_DOSAGE"]) <= args.vmiss_filter)
            )
            vmatch["PCA_ELIGIBLE"] = PCA_ELIGIBLE
            if PCA_ELIGIBLE is True:
                n_PCA_ELIGIBLE += 1

            if n_matched == 1:
                writer = csv.DictWriter(
                    csvfile, fieldnames=vmatch.keys(), delimiter="\t"
                )
                writer.writeheader()
            writer.writerow(vmatch)
    logger.info(
        "{}/{} ({:.2f}%) of TARGET variants matched the REFERENCE data".format(
            n_matched, count_var_t, 100 * n_matched / count_var_t
        )
    )
    logger.info(
        "{}/{} ({:.2f}%) of matched variants are eligible for PCA".format(
            n_PCA_ELIGIBLE, n_matched, 100 * n_PCA_ELIGIBLE / n_matched
        )
    )

    # Output counts
    logger.info(
        "Outputting variant counts -> intersect_counts_{}.txt".format(args.filter_chrom)
    )
    with open(
        outdir / "intersect_counts_{}.txt".format(args.filter_chrom), "w"
    ) as outf:
        outf.write("\n".join(map(str, [count_var_t, count_var_r, n_matched])))


def read_var_general(path, chrom=None):
    """
    General function for reading variant files from plink2 outputs
    :param path: path to variant file
    :param chrom: filter to specific chromosome
    :return: row of a df as a dict
    """
    with xopen(path, "rt") as f:
        if "bim" in pathlib.Path(path).name:
            reader = csv.reader(f, delimiter="\t")
            # yes, A1/A2 in bim isn't ref/alt
            fields = ["#CHROM", "ID", "pos_cm", "POS", "REF", "ALT"]
            if (chrom is None) or (chrom == "ALL"):
                for row in reader:
                    yield dict(zip(fields, row, strict=True))
            else:
                for row in reader:
                    row = dict(zip(fields, row, strict=True))
                    if row["#CHROM"] == chrom:
                        yield row
        else:
            fields = None
            if (chrom is None) or (chrom == "ALL"):
                for row in f:
                    if row.startswith("##"):
                        continue
                    else:
                        row = row.strip().split("\t")
                        if fields is None:
                            fields = row
                        else:
                            yield dict(zip(fields, row, strict=True))
            else:
                for row in f:
                    if row.startswith("##"):
                        continue
                    else:
                        row = row.strip().split("\t")
                        if fields is None:
                            fields = row
                        else:
                            row = dict(zip(fields, row, strict=True))
                            if row["#CHROM"] == chrom:
                                yield row


def sorted_join_variants(path_ref, path_target):
    f1_iter = read_var_general(path_ref)
    f2_iter = read_var_general(path_target)

    prev_key1 = None  # Initialize previous key for file 1
    prev_key2 = None  # Initialize previous key for file 2

    line1 = next(f1_iter, None)
    line2 = next(f2_iter, None)

    while line1 is not None and line2 is not None:
        key1 = line1["CHR:POS:A0:A1"]
        key2 = line2["CHR:POS:A0:A1"]

        # Check if lines are sorted by the key for each file
        if prev_key1 is not None and key1 < prev_key1:
            raise ValueError("REFERENCE keys are not sorted")
        if prev_key2 is not None and key2 < prev_key2:
            raise ValueError("TARGET keys are not sorted")

        if key1 == key2:
            line1.update(line2)
            yield line1
            prev_key1 = key1  # Update previous key for file 1
            prev_key2 = key2  # Update previous key for file 2

            line1 = next(f1_iter, None)
            line2 = next(f2_iter, None)
        elif key1 < key2:
            prev_key1 = key1  # Update previous key for file 1
            line1 = next(f1_iter, None)
        else:
            prev_key2 = key2  # Update previous key for file 2
            line2 = next(f2_iter, None)


def allele_complement(s):
    """
    Complement alleles
    :param s: allele to be complemented
    :return: complement
    """
    return (
        s.replace("A", "V")
        .replace("T", "X")
        .replace("C", "Y")
        .replace("G", "Z")
        .replace("V", "T")
        .replace("X", "A")
        .replace("Y", "G")
        .replace("Z", "C")
    )


def aaf2maf(aaf):
    """
    Convert alternative allele frequency (AAF) to minor allele frequency (MAF)
    :param aaf: alternative allele frequency
    :return: minor allele frequency (MAF)
    """
    if aaf > 0.5:
        return 1 - aaf
    else:
        return aaf


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description=_description_text(),
        epilog=_epilog_text(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-r",
        "--reference",
        dest="reference",
        required=True,
        help="path/to/REFERENCE/pvar",
    )
    parser.add_argument(
        "-t",
        "--target",
        dest="target",
        required=True,
        nargs="+",
        help="<Required> A list of paths of target genomic variants (.bim/pvar format). The .afreq and .vmiss files are also required for these files.",
    )
    parser.add_argument(
        "-c",
        "--chrom",
        dest="filter_chrom",
        required=False,
        help="Whether to limit matches to specific chromosome of the reference",
    )
    parser.add_argument(
        "--maf_target",
        dest="maf_filter",
        default=0,
        type=float,
        required=False,
        help="Filter: Minimum minor Allele Frequency for PCA eligibility",
    )
    parser.add_argument(
        "--geno_miss",
        dest="vmiss_filter",
        default=1,
        type=float,
        required=False,
        help="Filter: Maximum Genotype missingness for PCA eligibility",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="<Optional> Extra logging information",
    )
    parser.add_argument(
        "--outdir", dest="outdir", required=True, help="<Required> Output directory"
    )
    return parser.parse_args(args)


def _description_text() -> str:
    return textwrap.dedent(
        """\
    Program to find matched variants (same strand) between a set of reference and target data .pvar/bim files. This 
    also evaluate whether the variants in the TARGET are suitable for inclusion in  a PCA analysis (excludes strand 
    ambiguous and multi-allelic/INDEL variants), and can also uses the .afreq and .vmiss files exclude variants with 
    missingness and MAF filters. 
   """
    )


def _epilog_text() -> str:
    return textwrap.dedent("""""")


if __name__ == "__main__":
    run_intersect()
