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


    # Process target variants
    with xopen('target_variants.txt', 'wt') as outf:
        outf.write('CHR:POS:A0:A1\tID_TARGET\tREF_TARGET\tIS_MA_TARGET\tALT_FREQ\tF_MISS_DOSAGE\n')
        for path in args.target:
            pvar = read_var_general(path)

            loc_afreq = path.replace('.pvar.zst', '.afreq.gz')
            afreq = read_var_general(loc_afreq)

            loc_vmiss = path.replace('.pvar.zst', '.vmiss.gz')
            vmiss = read_var_general(loc_vmiss)

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


def read_var_general(path):
    with xopen(path, "rt") as f:
        # pvars do have a header column and support arbitrary columns
        reader = csv.DictReader(filter(lambda row: row[:2]!='##', f), delimiter="\t") # need to remove comments of VCF-like characters
        for row in reader:
            yield row
def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description=_description_text(),
        epilog=_epilog_text(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-t",
        "--target",
        dest="target",
        required=True,
        nargs="+",
        help="<Required> A list of paths of target genomic variants (.bim/pvar format)",
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
