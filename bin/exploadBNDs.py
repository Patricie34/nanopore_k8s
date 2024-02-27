#! /usr/bin/env python
import argparse
from cyvcf2 import VCF, Writer
import numpy as np
import re

BNDdict = {
    "[N": 1,
    "N]": 2,
    "]N": 3,
    "N[": 4,
}


def getBNDdict(breakpoint: str) -> str:
    # from breakpoint, e.g. ']7:152321850]N' return its type, e.g. 'DUP-like'
    splited = breakpoint.split(':')[0]
    try:
        replaced = re.sub(r'[a-zA-Z0-9]+', 'N', splited)
        return BNDdict[replaced[0:2]]
    except:
        return '?'


# # Parse command line
parser = argparse.ArgumentParser(
    description='Expload BND records')
parser.add_argument('-v', '--vcf', metavar='input.vcf.gz',
                    required=True, dest='vcf', help='input VCF file (required)')
parser.add_argument('-o', '--out', metavar='output.vcf.gz',
                    required=True, dest='outname', help='output VCF file (required)')


args = parser.parse_args()
vcf = VCF(args.vcf)
w = Writer(args.outname, vcf)

# vcf = VCF('/Volumes/share/share/110000-MED/110323-imgg/111323-plevova/01.NanoBreak/data/samples/BRNO3024/nano/VarCal/BRNO3024.svim.BNDs.vcf')
# w = Writer("TEST_INFLATE_DELLY.vcf", vcf)

patternChr = r'[\dA-Za-z]'
patternPos = r'[0-9]'
for record in vcf:
    w.write_record(record)  # original record
    chr1 = record.CHROM
    pos1 = record.POS
    chr2, pos2 = record.ALT[0].strip("NCTGA[]").split(":")
    record.ID = record.ID + 'dup'
    record.CHROM = chr2
    try:
        record.set_pos(int(pos2)-1)  # -1?
    except Exception as e:
        print(record.ALT)

    parts = record.ALT[0].split(':')
    BNDtype = getBNDdict(parts[0])
    if (BNDtype == 3):
        ALTnew = 'N[' + chr1 + ':' + str(pos1) + '['
    elif (BNDtype == 4):
        ALTnew = ']' + chr1 + ':' + str(pos1) + ']N'
    else:
        ALTnew = '[' + chr1 + ':' + str(pos1) + '[N'

    # match = re.search(patternChr, parts[0])
    # ALTchr = parts[0][:match.start()]+chr1

    # matches = list(re.finditer(patternPos, parts[1]))
    # ALTpos = str(pos1)+parts[1][matches[-1].end():]

    # ALTnew = ':'.join([ALTchr, ALTpos])
    record.ALT = [ALTnew]
    # print(record)
    w.write_record(record)  # duplicated record
