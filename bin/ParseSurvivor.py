import allel
import re
import csv

callset = allel.read_vcf('/Volumes/share/share/110000-MED/110323-imgg/111323-plevova/01.NanoBreak/data/samples/BRNO0627/nano/VarCal/SURVIVOR.Annotated.vcf',
                         fields=['samples', 'variants/ALT', 'variants/CHROM', 'variants/FILTER_PASS', 'variants/ID', 'variants/POS', 'variants/QUAL', 'variants/REF', 'QV', 'TY', 'CO', 'CSQ'])

sorted(callset.keys())
# callset['samples'][1]

# callset['calldata/GT'][0]

# callset['variants/ALT'][226]
# callset['variants/ID'][226]
# callset['calldata/QV'][226]
callset['variants/CSQ'][226].split('|')[4]

# callset['calldata/TY'][226] == 'TRA'

# if (callset['calldata/TY'][226].all() == 'TRA'):


# callset['calldata/TY']

for row in callset:
    print(row['calldata/TY'])


vcflen = len(callset['calldata/TY'])

BNDdict = {
    "[N": "INV-like",
    "N]": "INV-like",
    "N[": "DEL-like",
    "]N": "DUP-like"
}


def getBNDdict(breakpoint: str) -> str:
    # from breakpoint, e.g. ']7:152321850]N' return its type, e.g. 'DUP-like'
    splited = breakpoint.split(':')[0]
    replaced = re.sub(r'[a-zA-Z0-9]+', 'N', splited)
    return BNDdict[replaced[0:2]]


vars = []
for i in range(1, vcflen):
    if (callset['calldata/TY'][i].all() == 'TRA'):
        from_chr = str(callset['variants/CHROM'][i])
        from_pos = int(callset['variants/POS'][i])
        to_chr = str(callset['calldata/CO'][i][0].split('-')[1].split('_')[0])
        to_pos = int(callset['calldata/CO'][i][0].split('-')[1].split('_')[1])

        vars.append({
            "id": str(callset['variants/ID'][i]),
            "from_chr": from_chr,
            "from_pos": from_pos,
            "to_chr": to_chr,
            "to_pos": to_pos,
            "dist": "-" if (from_chr != to_chr) else str(abs(int(from_pos)-int(to_pos))),
            "BNDtype": getBNDdict(callset['variants/ALT'][i][0]),
            # "CNV_DB_ids": CNV_DB_ids,
            # "CNV_DB_count": CNV_DB_count,
            # "CNV_DB_types": CNV_DB_types,
            #  "reads": f"{rec}".split("READS=",1)[1].split("\t", 1)[0].split(","),
            #  "support" : len(f"{rec}".split("READS=",1)[1].split("\t", 1)[0].split(",")),
            "support": int(callset['calldata/QV'][i][0]),
            #   "coverage" : sum([x[0] for x in coverage]),
            # "percentage": len(reads)/int(coverage) if coverage != "?" else "?",
            # "coverage": coverage,
            # "cyto_from": CytoLoc,
            # "cyto_to": cyto_to,
            # "from_unique": f"{rec}".split("UNIQUE=")[1].split(";")[0],
            "from_variant": str(callset['variants/CSQ'][i].split('|')[1]),
            "from_gene": str(callset['variants/CSQ'][i].split('|')[3]),
            "from_gene_ID": str(callset['variants/CSQ'][i].split('|')[4]),
            # "to_unique": "?",
            # "to_variant": "?",
            # "to_gene": "?",
            # "to_gene_ID": "?",
            # "line": '\t'.join(f"{rec}".split('\t')[2:])
        })


csv_file_path = "output.csv"
# Extracting column names from the first dictionary in the list
fieldnames = vars[0].keys()
# Writing to CSV
with open(csv_file_path, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    # Write header
    writer.writeheader()

    # Write data
    writer.writerows(vars)
