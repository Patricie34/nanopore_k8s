import sys
from pysam import VariantFile

# VCFFile = "./SVIM_workdir/variants.vcf"

def read_variants(file_path: str) -> list:
  print(file_path)
  vcf = VariantFile(file_path[0])
  vars = []
  id = 0

  for rec in vcf.fetch():
    if "BND" in rec.id:
      [to_chr, to_pos] = rec.alleles[1].strip("N").strip("]").strip("[").split(":")
      from_chr = rec.chrom
      from_pos = rec.pos
      id += 1
      vars.append({ "id": id,
                    "from_chr": from_chr,
                     "from_pos": int(from_pos),
                     "to_chr": to_chr,
                     "to_pos": int(to_pos),
                     "dist": "-" if (from_chr != to_chr) else abs(int(from_pos)-int(to_pos)),
                     "match_snif": "",
                     "match_svim": "",
                     "match_cute": "",
                     "line": f"{rec}".strip()})

  return(vars)


def filter_duplicated(DupList):
  for i in range(len(DupList)-1):
    for j in range(i+1,len(DupList)-1):
      for read in DupList[i]["reads"]:
        if read in DupList[j]["reads"] and DupList[i]["from_pos"] == DupList[j]["to_pos"] and DupList[i]["from_chr"] == DupList[j]["to_chr"]:
          DupList[j]["reads"].remove(read)
          #if not len(DupList[j]["reads"]):
           # DedupList.pop(j)
  DedupList =  [rec for rec in DupList if len(rec["reads"])]
  return(DedupList)


if __name__ == "__main__":
   vars = read_variants(file_path = sys.argv[1:])
   filtered = filter_duplicated(vars)

   print("id\tdistance\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n")
   for v in filtered:
    print(f"{v['id']}\t{v['dist']}\t{v['line']}\n")