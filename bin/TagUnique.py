import sys
import pandas
from cyvcf2 import VCF, Writer


def tag_variants(file_path: str, uniques: str, name: str) -> list:
  print(file_path)
  ids = list(pandas.read_csv(uniques, header = None)[0])
  vcf = VCF(file_path)
  vcf.add_info_to_header({
    'ID': 'UNIQUE',
    'Description': 'Numeric code for filtering reason',
    'Type': 'String',
    'Number':'1'
})

 # create a new vcf Writer using the input vcf as a template.
  fname = f"{name}.UniqueTag.vcf"
  w = Writer(fname, vcf)

  for v in vcf:
     # The get_gene_intersections function is not shown.
     # This could be any operation to find intersections
     # or any manipulation required by the user.
   if v.ID in ids:
    v.INFO["UNIQUE"] = "FALSE"
   else:
    v.INFO["UNIQUE"] = "TRUE"
   w.write_record(v)

  w.close(); vcf.close()


if __name__ == "__main__":
   print("Taging unique vcf entries")
   tagged = tag_variants(file_path = sys.argv[1], uniques = sys.argv[2], name = sys.argv[3])