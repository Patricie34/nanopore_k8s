import sys
from pysam import VariantFile

# VCFFile = "./SVIM_workdir/variants.vcf"

def read_variants(file_path: str) -> list:
  vcf = VariantFile(file_path)
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


# with open('./SVIM_workdir/Variants_distances.tsv','w') as f:
#   f.write("id\tdistance\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n")
#   for v in vars:
#     f.write(f"{v['id']}\t{v['dist']}\t{v['line']}\n")


if __name__ == "__main__":
   vars = read_variants(file_path = sys.argv[1:])

   print("id\tdistance\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n")
   for v in vars:
    print(f"{v['id']}\t{v['dist']}\t{v['line']}\n")