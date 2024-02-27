#! /usr/bin/env python

import argparse
from cyvcf2 import VCF, Writer
import numpy as np


# Parse command line
parser = argparse.ArgumentParser(
    description='Filter singleton SVs from multi-sample Survivor .vcf')
parser.add_argument('-v', '--vcf', metavar='input.vcf.gz',
                    required=True, dest='vcf', help='input VCF file (required)')
parser.add_argument('-n', '--name', metavar='BRNOxxxx', required=True,
                    dest='sample', help='Sample name to be filtered out (required)')
parser.add_argument('-o', '--out', metavar='output.vcf.gz',
                    required=True, dest='outname', help='output VCF file (required)')


args = parser.parse_args()

# vcf = VCF('/Volumes/share/share/110000-MED/110323-imgg/111323-plevova/01.NanoBreak/data/samples/BRNO3024/nano/VarCal/survivor/merged.bcf')
# vcf.set_samples(['Sample', 'BRNO3024'])

vcf = VCF(args.vcf)
vcf.set_samples([args.sample+'.Sniffles', args.sample])

vcf.add_format_to_header({
    'ID': 'NUM_CALLERS',
    'Description': 'Number of variant callers supporting this variant',
    'Type': 'String',
    'Number': '1'
})

idx = vcf.samples.index(args.sample)
idxRef = vcf.samples.index(args.sample+'.Sniffles')
w = Writer(args.outname, vcf)
# idx = vcf.samples.index("BRNO3024")
# idxRef = vcf.samples.index('Sample')
# w = Writer('test.vcf', vcf)


for record in vcf:
    # Check if variant is present in filtered sample and is the only carrier
    # and (record.INFO.get('SUPP') == '1')
    if ((record.format("CO")[idx] != 'NAN')):
        # append information about how many variant callers detected this variant
        try:
            supp_value = record.format('NUM_CALLERS')[idxRef]
        except Exception as e:
            print(
                f"{record}")
            print(f"Caught an exception: {e}")
            supp_value = '?'
            # continue
        record.set_format('NUM_CALLERS', np.full(
            2, supp_value, dtype='S'))
        w.write_record(record)
