import allel
import pandas as pd
import numpy as np
import argparse


parser = argparse.ArgumentParser(
    description='Edit Tobias 1K ONT .tsv file')
parser.add_argument('-v', '--vcf', metavar='input.vcf',
                    required=True, dest='vcf', help='input VCF file (required)')
parser.add_argument('-t', '--tsv', metavar='classification from Tobias',
                    required=True, dest='tsv', help='classification from Tobias')
args = parser.parse_args()


def flatten_nested_dict(d):
    flattened_dict = {}
    for key, value in d.items():
        # flatten nested array, take first item in the nested array
        if isinstance(value[0], np.ndarray):
            flattened_dict[key] = [sublist[0] for sublist in value]
        else:
            flattened_dict[key] = value
    return flattened_dict


vcf = allel.read_vcf(args.vcf, fields=["ID", "ALT"])
flattened = flatten_nested_dict(vcf)
callset_pd = pd.DataFrame(flattened)

tsv_pd = pd.read_csv(args.tsv, delimiter='\t')

merged_df = pd.merge(tsv_pd, callset_pd,
                     right_on='variants/ID',
                     left_on="ID",
                     how='left')
print("id", "chr", "pos", "chr2", "pos2",
      "svtype", "svlen", "CT", sep="\t")

for i, row in merged_df.iterrows():
    chr_to = row['variants/ALT'].strip("NCTGA[]").split(":")[0] if row['ID'].startswith("BND") else '-'
    print(row.ID, row.Chrom, row.Start, chr_to, row.End,
          row.SVType, row.Length, row.CT, sep='\t')
