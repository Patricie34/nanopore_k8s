import sys
import pysam

#file_path = "/Volumes/lamb/shared/MedGen/nanobreak/src/pipeline/project/xsvato01/nanopore_k8s/bin/500.Unique.Annotated.vcf"
#bam_path = "/Volumes/share/share/110000-MED/999999-gic/01.NanoBreak/data/samples/BRNO0627/nano/mapped/BRNO0627.sorted.bam"

def read_variants(file_path: str, bam_path: str) -> list:
  print(file_path)
  vcf = pysam.VariantFile(file_path)
  bam = pysam.AlignmentFile(bam_path, mode = 'rb')  # 'rb' ~ read bam

  vars = []
  id = 0

  for rec in vcf.fetch():
    if "BND" in rec.id:
      [to_chr, to_pos] = rec.alleles[1].strip("N").strip("]").strip("[").split(":")
      from_chr = rec.chrom
      from_pos = rec.pos
      id += 1
      try:
        coverage = bam.count_coverage(
                    contig = from_chr,     # Chromosome ID; also might be "chr1" or similar 
                    start = from_pos-1,
                    stop = from_pos,
                    )
      except ValueError:
        print(f'Something went wrong at {rec}, from chr: {from_chr} from pos: {from_pos}.')
      
      vars.append({ "id": id,
                    "from_chr": from_chr,
                     "from_pos": int(from_pos),
                     "to_chr": to_chr,
                     "to_pos": int(to_pos),
                     "dist": "-" if (from_chr != to_chr) else str(abs(int(from_pos)-int(to_pos))),
                     "reads": f"{rec}".split("READS=",1)[1].split("\t", 1)[0].split(";")[0].split(","),
                     #"reads": f"{rec}".split("READS=",1)[1].split("\t", 1)[0].split(","),
                     #"support" : len(f"{rec}".split("READS=",1)[1].split("\t", 1)[0].split(",")),
                     "support" : len(f"{rec}".split("READS=",1)[1].split("\t", 1)[0].split(";")[0].split(",")),
                     "coverage" : sum([x[0] for x in coverage]),
                     "from_unique" : f"{rec}".split("UNIQUE=")[1].split(";")[0],
                     "from_variant" : f"{rec}".split('|')[1],
                     "from_gene": f"{rec}".split('|')[3],
                     "from_gene_ID" : f"{rec}".split('|')[4],
                     "to_unique" : "?",
                     "to_variant" : "?",
                     "to_gene" : "?",
                     "to_gene_ID" : "?",
                     "line": '\t'.join(f"{rec}".split('\t')[2:])})

  return(vars)


def filter_duplicated(DupList):
  listlen = len(DupList)-1
  for i in range(listlen):
    for j in range(i+1,listlen):
      for read in DupList[i]["reads"]:
        if read in DupList[j]["reads"] and DupList[i]["from_pos"] == DupList[j]["to_pos"] and DupList[i]["from_chr"] == DupList[j]["to_chr"]:
          DupList[j]["reads"].remove(read)
          DupList[i]["to_unique"] = DupList[j]['line'].split("UNIQUE=")[1].split(";")[0]
          DupList[i]["to_variant"] = DupList[j]['line'].split('|')[1]
          DupList[i]["to_gene"] = DupList[j]['line'].split('|')[3]
          DupList[i]["to_gene_ID"] = DupList[j]['line'].split('|')[4]

          #if not len(DupList[j]["reads"]):
            #DedupList.pop(j)
            
  DedupList =  [rec for rec in DupList if len(rec["reads"])]
  return(DedupList)


if __name__ == "__main__": #arg[1] = cesta k csv arg[2] = cesta k bam arg[3] = jmeno k ulozeni
   print("Parsing vcf")
   vars = read_variants(file_path = sys.argv[1], bam_path = sys.argv[2] )
   #vars = read_variants(file_path, bam_path )

   print("Filtering vcf")
   filtered = filter_duplicated(vars)

  # print("id\tdistance\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample")
   #for v in filtered:
    #print(f"{v['id']}\t{v['dist']}\t{v['line']}")
   print("Saving vcf")
   with open(f"Dedup_{sys.argv[3]}.tsv",'w') as f:
    f.write("id\tdistance\t"\
    "ChrFrom\tPosFrom\tFromUnique?\tFromVariant\tFromGene\tFromGeneID\t"\
    "ChrTo\tPosTo\tToUnique?\tToVariant\tToGene\tToGeneID\t"\
    "Support_nonUnique\tCoverage\t\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for v in filtered:
      f.write(\
      f"{v['id']}\t{v['dist']}\t"\
      f"{v['from_chr']}\t{v['from_pos']}\t{v['from_unique']}\t{v['from_variant']}\t{v['from_gene']}\t{v['from_gene_ID']}\t"\
      f"{v['to_chr']}\t{v['to_pos']}\t{v['to_unique']}\t{v['to_variant']}\t{v['to_gene']}\t{v['to_gene_ID']}\t"\
      f"{v['support']}\t{v['coverage']}\t\t{v['line']}")

   with open(f"Dedup_1000dist{sys.argv[3]}.tsv",'w') as f:
    f.write("id\tdistance\t"\
    "ChrFrom\tPosFrom\tFromUnique?\tFromVariant\tFromGene\tFromGeneID\t"\
    "ChrTo\tPosTo\tToUnique?\tToVariant\tToGene\tToGeneID\t"\
    "Support_nonUnique\tCoverage\t\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for v in filtered:
     if not v['dist'].isnumeric() or int(v['dist']) > 1000:
      f.write(\
      f"{v['id']}\t{v['dist']}\t"\
      f"{v['from_chr']}\t{v['from_pos']}\t{v['from_unique']}\t{v['from_variant']}\t{v['from_gene']}\t{v['from_gene_ID']}\t"\
      f"{v['to_chr']}\t{v['to_pos']}\t{v['to_unique']}\t{v['to_variant']}\t{v['to_gene']}\t{v['to_gene_ID']}\t"\
      f"{v['support']}\t{v['coverage']}\t\t{v['line']}")
