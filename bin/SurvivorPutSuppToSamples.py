from cyvcf2 import VCF, Writer
import numpy as np
import argparse


# Parse command line
parser = argparse.ArgumentParser(
    description='Add SUPP info to Samples so that after SURVIVOR private filtering the information about how many variant callers detected given variant stays')
parser.add_argument('-v', '--vcf', metavar='input.vcf',
                    required=True, dest='vcf', help='input VCF file (required)')
parser.add_argument('-o', '--outname', metavar='output.vcf',
                    required=True, dest='out', help='output VCF file (required)')

args = parser.parse_args()

# Parse VCF

vcf = VCF(args.vcf)
# create a new vcf Writer using the input vcf as a template.
fname = "out.vcf"

vcf.add_format_to_header({
    'ID': 'NUM_CALLERS',
    'Description': 'Number of variant callers supporting this variant',
    'Type': 'String',
    'Number': '1'
})

w = Writer(args.out, vcf)
nsamples = len(vcf.samples)
for v in vcf:
    # Retrieve the SUPP value from the INFO field
    supp_value = v.INFO.get('SUPP')
    v.set_format('AAL', np.full(
        nsamples, supp_value, dtype='S'))
    v.set_format('NUM_CALLERS', np.full(
        nsamples, supp_value, dtype='S'))
    w.write_record(v)


w.close()
vcf.close()
