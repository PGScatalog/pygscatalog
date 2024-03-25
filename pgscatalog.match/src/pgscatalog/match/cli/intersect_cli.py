import argparse
import logging
import sys
from xopen import xopen
import csv
import textwrap


from pgscatalog.core import TargetVariants

logger = logging.getLogger(__name__)


def run_intersect():
    args = parse_args()

    # Process reference variants
    with xopen('reference_variants.txt', 'wt') as outf:
        outf.write('CHR:POS:A0:A1\tID_REF\tREF_REF\tIS_INDEL\tSTRANDAMB\tIS_MA_REF\n')
        ref_pvar = read_var_general(args.reference, chrom=args.filter_chrom)
        for v in ref_pvar:
            ALTs = v['ALT'].split(',')
            IS_MA_REF = len(ALTs) > 1
            for i, ALT in enumerate(ALTs):
                if v['REF'] < ALT:
                    key = '{}:{}:{}:{}'.format(v['#CHROM'], v['POS'], v['REF'], ALT)
                else:
                    key = '{}:{}:{}:{}'.format(v['#CHROM'], v['POS'], ALT, v['REF'])

                IS_INDEL = (len(v['REF']) > 1) | (len(ALT) > 1)
                STRANDAMB = (v['REF'] == allele_complement(ALT))
                outf.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(key,v['ID'], v['REF'], IS_INDEL, STRANDAMB, IS_MA_REF))

    # Process target variants
    with xopen('target_variants.txt', 'wt') as outf:
        outf.write('CHR:POS:A0:A1\tID_TARGET\tREF_TARGET\tIS_MA_TARGET\tALT_FREQ\tF_MISS_DOSAGE\n')
        for path in args.target:
            pvar = read_var_general(path)

            loc_afreq = path.replace('.pvar.zst', '.afreq.gz')
            afreq = read_var_general(loc_afreq)

            loc_vmiss = path.replace('.pvar.zst', '.vmiss.gz')
            vmiss = read_var_general(loc_vmiss)

            ## ToDo: figure out if it needs to change to bim files
            for v, freq, miss in zip(pvar, afreq, vmiss):
                # if v['ID'] != freq['ID'] != miss['ID']:
                #     print(v)
                ALTs = v['ALT'].split(',')
                ALT_FREQS = freq['ALT_FREQS'].split(',')
                F_MISS_DOSAGE = miss['F_MISS_DOSAGE']
                IS_MA_TARGET = len(ALTs) > 1
                for i, ALT in enumerate(ALTs):
                    if v['REF'] < ALT:
                        key = '{}:{}:{}:{}'.format(v['#CHROM'], v['POS'], v['REF'], ALT)
                    else:
                        key = '{}:{}:{}:{}'.format(v['#CHROM'], v['POS'], ALT, v['REF'])
                    outf.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(key,v['ID'],v['REF'], str(IS_MA_TARGET), ALT_FREQS[i], F_MISS_DOSAGE))


def read_var_general(path, chrom=None):
    with xopen(path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            else:
                fieldnames = line.strip().split("\t")
        reader = csv.DictReader(f, fieldnames=fieldnames, delimiter="\t")
        if chrom is None:
            for row in reader:
                yield row
        else:
            for row in reader:
                if row['#CHROM'] == chrom:
                    yield row


def allele_complement(s):
    return s.replace("A", "V").replace("T", "X").replace("C", "Y").replace("G", "Z").replace("V", "T").replace("X", "A").replace("Y", "G").replace("Z", "C")

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
        help="<Required> A list of paths of target genomic variants (.bim/pvar format)",
    )
    parser.add_argument(
        "-c",
        "--chrom",
        dest="filter_chrom",
        required=False,
        help="whether to limit matches to specific chromosome of the reference",
    )
    return parser.parse_args(args)


def _description_text() -> str:
    return textwrap.dedent(
        """\
    PLACEHOLDER
   """
    )


def _epilog_text() -> str:
    return textwrap.dedent(
        """\
    PLACEHOLDER
    """
    )
