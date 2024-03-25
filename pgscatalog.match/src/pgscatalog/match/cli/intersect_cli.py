import argparse
import logging
from xopen import xopen
import csv
import textwrap
import heapq


logger = logging.getLogger(__name__)


def run_intersect():
    args = parse_args()

    # Process & sort reference variants
    with xopen('reference_variants.txt', 'wt') as outf:
        outf.write('CHR:POS:A0:A1\tID_REF\tREF_REF\tIS_INDEL\tSTRANDAMB\tIS_MA_REF\n')
        ref_heap = []
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
                heapq.heappush(ref_heap, ([key, v['ID'], v['REF']],[IS_INDEL, STRANDAMB, IS_MA_REF]))

        # Output the sorted reference variants
        n_ref = len(ref_heap)
        for i in range(n_ref):
            popped = heapq.heappop(ref_heap)
            outf.write('\t'.join([str(x) for x in popped[0] + popped[1]]) + '\n')
        del ref_heap

    # Process & sort target variants
    # ToDo: check if it works for bim format files?
    with xopen('target_variants.txt', 'wt') as outf:
        outf.write('CHR:POS:A0:A1\tID_TARGET\tREF_TARGET\tIS_MA_TARGET\tALT_FREQ\tF_MISS_DOSAGE\n')
        target_heap = []
        for path in args.target:
            pvar = read_var_general(path, chrom=None)  # essential not to filter if it is target (messes up common line indexing)

            loc_afreq = path.replace('.pvar.zst', '.afreq.gz')
            afreq = read_var_general(loc_afreq, chrom=None)  # essential not to filter if it is target (messes up common line indexing)

            loc_vmiss = path.replace('.pvar.zst', '.vmiss.gz')
            vmiss = read_var_general(loc_vmiss, chrom=None)  # essential not to filter if it is target (messes up common line indexing)

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
                    # outf.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(key, v['ID'], v['REF'], str(IS_MA_TARGET), ALT_FREQS[i],
                    #                                              F_MISS_DOSAGE))
                    heapq.heappush(target_heap, ([key, v['ID'], v['REF']], [IS_MA_TARGET, ALT_FREQS[i],F_MISS_DOSAGE]))

        # Output the sorted reference variants
        n_target = len(target_heap)
        for i in range(n_target):
            popped = heapq.heappop(target_heap)
            outf.write('\t'.join([str(x) for x in popped[0] + popped[1]]) + '\n')
        del target_heap

    # Merge matched variants on sorted files
    n_matched = 0
    with open('matched_variants.txt', 'w') as csvfile:
        for vmatch in sorted_join_variants('reference_variants.txt', 'target_variants.txt'):
            n_matched += 1
            vmatch['SAME_REF'] = vmatch['REF_REF'] == vmatch['REF_REF']
            if n_matched == 1:
                writer = csv.DictWriter(csvfile, fieldnames=vmatch.keys(), delimiter='\t')
                writer.writeheader()
            writer.writerow(vmatch)

    # Output counts
    with open('intersect_counts_${}.txt'.format(chrom), 'w') as outf:
        outf.write('\n'.join(map(str, [n_target, n_ref, n_matched])))


def read_var_general(path, chrom=None):
    with xopen(path, "rt") as f:
        # ToDo: check if this is memory inefficent
        reader = csv.DictReader(filter(lambda row: row[:2]!='##', f), delimiter="\t") # need to remove comments of VCF-like characters, might be fully in memory though
        if (chrom == None) or (chrom == 'ALL'):
            for row in reader:
                yield row
        else:
            for row in reader:
                if row['#CHROM'] == chrom:
                    yield row


def sorted_join_variants(reffile, targetfile):
    f1_iter = read_var_general(reffile)
    f2_iter = read_var_general(targetfile)

    line1 = next(f1_iter, None)
    line2 = next(f2_iter, None)

    while line1 is not None and line2 is not None:
        key1 = line1['CHR:POS:A0:A1']
        key2 = line2['CHR:POS:A0:A1']

        if key1 == key2:
            line1.update(line2)
            yield line1
            line1 = next(f1_iter, None)
            line2 = next(f2_iter, None)
        elif key1 < key2:
            line1 = next(f1_iter, None)
        else:
            line2 = next(f2_iter, None)


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
