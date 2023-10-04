import sys
import pysam
import time
import pandas
import re


BNDdict = {
  "[N": "INV-like",
  "N]": "INV-like",
  "N[": "DEL-like",
  "]N": "DUP-like"
}

def getBNDdict(breakpoint: str) -> str:
    #from breakpoint, e.g. ']7:152321850]N' return its type, e.g. 'DUP-like'
    splited = breakpoint.split(':')[0]
    replaced = re.sub(r'[a-zA-Z0-9]+', 'N', splited)
    return BNDdict[replaced[0:2]]

#file_path = "/Volumes/share/share/110000-MED/110323-imgg/111323-plevova/01.NanoBreak/data/samples/BRNO2013/nano/mapped/Annotated5000.vcf"
#bam_path = "/Volumes/share/share/110000-MED/999999-gic/01.NanoBreak/data/samples/BRNO0627/nano/mapped/BRNO0627.sorted.bam"
#coverage_path = "/Volumes/share/share/110000-MED/110323-imgg/111323-plevova/01.NanoBreak/data/samples/BRNO2013/nano/mapped/Coverage.txt"


def read_variants(file_path: str, coverage_path: str) -> list:
    vcf = pysam.VariantFile(file_path)
    bam_cov = pandas.read_csv(
        coverage_path, delimiter="\t", header=None, index_col=False)
    bam_cov.columns = ['Chr', 'Start', 'Stop', 'ID', 'Chujovina', 'Coverage', 'CytoLoc','CNV_DB_ids','CNV_DB_count','CNV_DB_types']
    vars = []
    id = 0

    for rec in vcf.fetch():
        if "BND" in rec.id:
            [to_chr, to_pos] = rec.alleles[1].strip(
                "N").strip("]").strip("[").split(":")
            from_chr = rec.chrom
            from_pos = rec.pos
            CytoLoc = bam_cov[bam_cov["ID"] ==
                                   rec.id]["CytoLoc"].to_string(index=False)
            CNV_DB_ids = bam_cov[bam_cov["ID"] ==
                                   rec.id]["CNV_DB_ids"].to_string(index=False)
            CNV_DB_count = bam_cov[bam_cov["ID"] ==
                                   rec.id]["CNV_DB_count"].to_string(index=False)
            CNV_DB_types = bam_cov[bam_cov["ID"] ==
                                   rec.id]["CNV_DB_types"].to_string(index=False)
            id += 1
            if from_chr == "MT" or to_chr == "MT":
                print(f"Skipping {rec.id}")
                continue
            try:
                coverage = bam_cov[bam_cov["ID"] ==
                                   rec.id]["Coverage"].to_string(index=False)
            except ValueError:
                print(
                    f"Something went wrong at {rec}, from chr: {from_chr} from pos: {from_pos}")
                coverage = "?"
                reads = f"{rec}".split("READS=", 1)[1].split("\t", 1)[0].split(";")[0].split(",")

            vars.append({"id": id,
                         "from_chr": from_chr,
                         "from_pos": int(from_pos),
                         "to_chr": to_chr,
                         "to_pos": int(to_pos),
                         "dist": "-" if (from_chr != to_chr) else str(abs(int(from_pos)-int(to_pos))),
                         "reads": reads,
                         "BNDtype": getBNDdict(rec.alleles[1]),
                         "CNV_DB_ids": CNV_DB_ids,
                         "CNV_DB_count": CNV_DB_count,
                         "CNV_DB_types": CNV_DB_types,
                        #  "reads": f"{rec}".split("READS=",1)[1].split("\t", 1)[0].split(","),
                         #  "support" : len(f"{rec}".split("READS=",1)[1].split("\t", 1)[0].split(",")),
                         "support": len(reads),
                         #   "coverage" : sum([x[0] for x in coverage]),
                         "percentage": len(reads)/coverage,
                         "coverage": coverage,
                         "cyto_map" : CytoLoc,
                         "from_unique": f"{rec}".split("UNIQUE=")[1].split(";")[0],
                         "from_variant": f"{rec}".split('|')[1],
                         "from_gene": f"{rec}".split('|')[3],
                         "from_gene_ID": f"{rec}".split('|')[4],
                         "to_unique": "?",
                         "to_variant": "?",
                         "to_gene": "?",
                         "to_gene_ID": "?",
                         "line": '\t'.join(f"{rec}".split('\t')[2:])})

    return (vars)

def filter_duplicated(DupList):
    listlen = len(DupList)-1
    for i in range(listlen):
        for j in range(i+1, listlen):
         # for read in DupList[i]["reads"]:
            #        if read in DupList[j]["reads"] and DupList[i]["from_pos"] == DupList[j]["to_pos"] and DupList[i]["from_chr"] == DupList[j]["to_chr"]:
         #       print(j)
            if DupList[i]["from_pos"] == DupList[j]["to_pos"] and DupList[i]["from_chr"] == DupList[j]["to_chr"]:

                # DupList[j]["reads"].remove(read)
                DupList[j]["reads"] = False
                DupList[i]["to_unique"] = DupList[j]['line'].split("UNIQUE=")[
                    1].split(";")[0]
                DupList[i]["to_variant"] = DupList[j]['line'].split('|')[1]
                DupList[i]["to_gene"] = DupList[j]['line'].split('|')[3]
                DupList[i]["to_gene_ID"] = DupList[j]['line'].split('|')[4]
                break

    # DedupList =  [rec for rec in DupList if len(rec["reads"])]
    DedupList = [rec for rec in DupList if (rec["reads"])]
    return (DedupList)


if __name__ == "__main__":
    # arg[1] = cesta k csv arg[2] = cesta k bam arg[3] = jmeno k ulozeni
    print("Parsing vcf")
    t = time.time()
    vars = read_variants(file_path=sys.argv[1], coverage_path=sys.argv[2])
    #vars = read_variants(file_path, coverage_path)
    print(time.time() - t)

    print("Filtering vcf")
    filtered = filter_duplicated(vars)
    print(time.time() - t)

    # print("id\tdistance\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample")
    # for v in filtered:
    # print(f"{v['id']}\t{v['dist']}\t{v['line']}")
    print("Saving vcf")
    with open(f"Dedup_{sys.argv[3]}.tsv", 'w') as f:
        f.write("id\tdistance\t"
                "ChrFrom\tPosFrom\tCytoMap\tBNDtype\tCNV_DB_ids\tCNV_DB_count\tCNV_DB_types\t"
                "FromUnique?\tFromVariant\tFromGene\tFromGeneID\t"
                "ChrTo\tPosTo\tToUnique?\tPercentage\tToVariant\tToGene\tToGeneID\t"
                "Support_nonUnique\tCoverage\t\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in filtered:
            f.write(
                f"{v['id']}\t{v['dist']}\t"
                f"{v['from_chr']}\t{v['from_pos']}\t{v['cyto_map']}\t{v['BNDtype']}\t{v['CNV_DB_ids']}\t{v['CNV_DB_count']}\t{v['CNV_DB_types']}\t"
                f"{v['from_unique']}\t{v['from_variant']}\t{v['from_gene']}\t{v['from_gene_ID']}\t"
                f"{v['to_chr']}\t{v['to_pos']}\t{v['to_unique']}\t{v['percentage']}\t{v['to_variant']}\t{v['to_gene']}\t{v['to_gene_ID']}\t"
                f"{v['support']}\t{v['coverage']}\t\t{v['line']}")

    with open(f"Dedup_1000dist{sys.argv[3]}.tsv", 'w') as f:
        f.write("id\tdistance\t"
                "ChrFrom\tPosFrom\tCytoMap\tBNDtype\tCNV_DB_ids\tCNV_DB_count\tCNV_DB_types\t"
                "FromUnique?\tFromVariant\tFromGene\tFromGeneID\t"
                "ChrTo\tPosTo\tToUnique?\tPercentage\tToVariant\tToGene\tToGeneID\t"
                "Support_nonUnique\tCoverage\t\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in filtered:
            if not v['dist'].isnumeric() or int(v['dist']) > 1000:
                f.write(
                    f"{v['id']}\t{v['dist']}\t"
                    f"{v['from_chr']}\t{v['from_pos']}\t{v['cyto_map']}\t{v['BNDtype']}\t{v['CNV_DB_ids']}\t{v['CNV_DB_count']}\t{v['CNV_DB_types']}\t"
                    f"{v['from_unique']}\t{v['from_variant']}\t{v['from_gene']}\t{v['from_gene_ID']}\t"
                    f"{v['to_chr']}\t{v['to_pos']}\t{v['to_unique']}\t{v['percentage']}\t{v['to_variant']}\t{v['to_gene']}\t{v['to_gene_ID']}\t"
                    f"{v['support']}\t{v['coverage']}\t\t{v['line']}")