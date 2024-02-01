#! /usr/bin/env python

from __future__ import print_function
import argparse
import sys
import cyvcf2


# Parse command line
parser = argparse.ArgumentParser(description='Filter singleton SVs from multi-sample BCF')
parser.add_argument('-v', '--vcf', metavar='input.vcf.gz', required=True, dest='vcf', help='input VCF file (required)')
args = parser.parse_args()

# Parse VCF
vcf = cyvcf2.VCF(args.vcf)
samples = list(vcf.samples)
print("chr", "pos", "chr2", "pos2", "id", "svtype", "svlen", "strand", "quality", "support", "carrier", sep="\t")
for record in vcf:
    # Ignore multi-allelics
    if len(record.ALT) > 1:
        continue

    # Random chromosome
    if record.CHROM.startswith("chrUn"):
        continue
    if record.CHROM.endswith("_random"):
        continue

    # One carrier?
    carrier = set()
    for s in samples:
        idx = samples.index(s)
        if ((record.gt_types[idx] != vcf.HOM_REF) and (record.gt_types[idx] != vcf.UNKNOWN)):
            carrier.add(s)
    if len(carrier) == 1:
        svtype = record.INFO.get("SVTYPE")
        support = record.INFO.get("SR")
        if support is None:
            support = record.INFO.get("SUPPORT")
        strand = record.INFO.get("CT")
        if strand is None:
            strand = record.INFO.get("STRAND")
        endval = record.INFO.get("END")
        svlen = record.INFO.get("SVLEN")
        if (svlen is None) and (endval is not None):
            svlen = endval - record.POS
        if (svlen is not None) and (svlen < 0):
            svlen = -1 * svlen
        if (svtype == "BND"):
            svlen = 0
            alt = record.ALT[0]
            alt = alt.replace('[','').replace(']','').replace('N','').replace('A','').replace('C','').replace('G','').replace('T','').split(':')
            chr2 = alt[0]
            pos2 = alt[1]
        else:
            chr2 = record.CHROM
            pos2 = endval
        if chr2.startswith("chrUn"):
            continue
        if chr2.endswith("_random"):
            continue
        print(record.CHROM, record.POS, chr2, pos2, record.ID, svtype, svlen, strand, record.QUAL, support, carrier.pop(), sep='\t')
        
