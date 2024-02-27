import allel
import re
import pandas as pd
import numpy as np
import argparse


def flatten_nested_dict(d):
    flattened_dict = {}
    n = len(d['calldata/NUM_CALLERS'])
    for key, value in d.items():
        # flatten nested array, take first item in the nested array
        if isinstance(value[0], np.ndarray):
            flattened_dict[key] = [sublist[0] for sublist in value]
        # prolong single items, eg sample name
        elif len(value) == 1:
            flattened_dict[key] = np.repeat(
                value, n)
        else:
            flattened_dict[key] = value

    return flattened_dict


parser = argparse.ArgumentParser(
    description='Expload BND records')
parser.add_argument('-v', '--vcf', metavar='input.vcf.gz',
                    required=True, dest='vcf', help='input VCF file (required)')
parser.add_argument('-c', '--cov', metavar='input.txt',
                    required=True, dest='cov', help='input cov file (required)')
parser.add_argument('-o', '--out', metavar='output.vcf.gz',
                    required=True, dest='outname', help='output VCF file (required)')


args = parser.parse_args()

callset_vcf = allel.read_vcf(args.vcf,
                             fields=['variants/CHROM', 'variants/POS', 'variants/SUPP', 'variants/SVLEN', 'calldata/NUM_CALLERS', 'calldata/CO', 'variants/ALT', 'variants/ID'])

# callset_vcf = allel.read_vcf(
#     '/Volumes/share/share/110000-MED/110323-imgg/111323-plevova/01.NanoBreak/data/samples/BRNO3024/nano/VarCal/survivor/singletons.vcf', fields=['*'])


flattened = flatten_nested_dict(callset_vcf)

vcf_pd = pd.DataFrame(flattened)


BNDdict = {
    "[N": "INV-like",
    "N]": "INV-like",
    "N[": "DEL-like",
    "]N": "DUP-like"
}


def getBNDdict(breakpoint: str) -> str:
    # from breakpoint, e.g. ']7:152321850]N' return its type, e.g. 'DUP-like'
    splited = breakpoint.split(':')[0]
    try:
        replaced = re.sub(r'[a-zA-Z0-9]+', 'N', splited)
        return BNDdict[replaced[0:2]]
    except:
        return '?'


dtype_cov = {
    'Chr': object,
    'Start': int,
    'Stop': int,
    'ID': object,
    'ALT': object,
    'Chujovina': object,
    'Coverage': int,
    'CytoLoc': object,
    'CNV_DB_ids': object,
    'CNV_DB_count': int,
    'CNV_DB_types': object
}

# Specify column names
column_names_cov = ['Chr', 'Start', 'Stop', 'ID', 'ALT', 'Chujovina',
                    'Coverage', 'CytoLoc', 'CNV_DB_ids', 'CNV_DB_count', 'CNV_DB_types']

# Read the CSV file with specified column names and data types
cov_file = pd.read_csv(args.cov,
                       delimiter="\t", header=None, names=column_names_cov, dtype=dtype_cov)

merged_df = pd.merge(vcf_pd, cov_file,
                     # left_on='variants/ID',
                     #  right_on="ID",
                     left_on=['variants/CHROM', 'variants/POS'],
                     right_on=['Chr', 'Stop'],
                     how='left')

vars_list = []

for i in range(1, len(vcf_pd)):
    try:
        splited_pos = merged_df['calldata/CO'][i].split('-')[1].split('_')
        # if (len(splited_pos) == 1):
        #     to_chr =
        # parsed_to = merged_df['calldata/CO'][i].split('-')[1].split('_')

        vars_list.append({
            "id": str(merged_df.at[i, 'variants/ID']),
            "from_chr": str(merged_df.at[i, 'variants/CHROM']),
            "from_pos": int(merged_df.at[i, 'variants/POS']),
            "to_chr": str(splited_pos[0]),
            "to_pos": int(splited_pos[1]),
            "cyto": merged_df.at[i, 'CytoLoc'],
            "dist": merged_df.at[i, 'variants/SVLEN'],
            "BNDtype": getBNDdict(merged_df.at[i, 'variants/ALT']),
            "CNV_DB_ids": merged_df.at[i, "CNV_DB_ids"],
            "CNV_DB_count": merged_df.at[i, 'CNV_DB_count'],
            "CNV_DB_types": merged_df.at[i, 'CNV_DB_types'],
            "SampleCount": int(merged_df.at[i, 'variants/SUPP']),
            "CallersCount": merged_df.at[i, 'calldata/NUM_CALLERS'],
            "coverage": int(merged_df.at[i, 'Coverage']),
        })
    except Exception as e:
        print(f"Caught an exception: {e}, {i}: {merged_df['calldata/CO'][i]}")
        continue

pd.DataFrame(vars_list).to_csv(args.out, sep='\t', index=False)
