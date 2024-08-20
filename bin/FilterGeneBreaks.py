import pandas as pd
import argparse


parser = argparse.ArgumentParser(
    description='Filter GeneBreaks.tsv from 1000g.tsv')
parser.add_argument('--geneBreaks', metavar='geneBreaks.tsv',
                    required=True, dest='geneBreaks', help='input geneBreaks file (required)')
parser.add_argument('--filteredGenomes', metavar='Input 1000g.tsv',
                    required=True, dest='filteredGenomes', help='filteredGenomes from Tobias (required)')
args = parser.parse_args()
geneBreaks = pd.read_csv(args.geneBreaks, delimiter='\t')

filteredGenomes = pd.read_csv(args.filteredGenomes, delimiter='\t')

finalMerged = pd.merge(filteredGenomes, geneBreaks,
                       right_on="query.id",
                       left_on="id",
                       how='left')

print("id", "chr", "pos", "chr2", "pos2" ,"refCov", "varCov",
      "svtype", "svlen", "CT", "quality", "startGene", "endGene", sep="\t")
for i, row in finalMerged.iterrows():
    print(row.id, row.chr, row.pos, row.chr2, row.pos2, row.refCov, row.varCov,
          row.svtype, row.svlen, row.CT, row['query.qual'], row['query.startfeature'], row['query.endfeature'], sep='\t')
