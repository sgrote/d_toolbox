#!/usr/bin/python3

''' 
subsample streamed in vcf-like file randomly given the proprtion of passing sites

output to stdout
input stream: (#header), columns: chr, pos, ...
input prob: proportion of passing sites; each position of the vcf will pass with the probability given by prob

'''


import sys
import argparse
from random import random


def main():
	parser = argparse.ArgumentParser(description="subsample streamed in vcf-like file randomly given the proprtion of passing sites.", usage="zcat file.vcf.gz | subsample_vcf.py 0.5")
	parser.add_argument("prob", type=float, help="proportion of passing sites; each position of the vcf will pass with the probability given by prob.")

	args = parser.parse_args()
	
	if sys.stdin.isatty():
		sys.exit("Error: Needs vcf from stdin, e.g. 'zcat file.vcf.gz | subsample_vcf.py 0.5'")
	vcf_file = sys.stdin
	
	# check prob
	if args.prob < 0 or args.prob > 1:
		sys.exit("prob has to be float in [0,1]")

	## print header
	vcf = vcf_file.readline()
	while vcf[0] == "#":
		sys.stdout.write(vcf)
		vcf = vcf_file.readline()
	
	## print position given pass-probability
	while not len(vcf) == 0:
		rannum = random() # random float [0,1)
		if rannum < args.prob:
			 sys.stdout.write(vcf)
		vcf = vcf_file.readline()


if __name__ == "__main__":
	main()	


''' test
SUBS=/mnt/expressions/steffi/D/d_toolbox/subsample_vcf.py
VCF=/mnt/scratch/steffi/D/Vcfs/mergedArchModernApes/merged_high_chr21.vcf.gz
zcat $VCF | wc -l
zcat $VCF | $SUBS 0.5 | wc -l
zcat $VCF | $SUBS 0.75 | wc -l
zcat $VCF | $SUBS 0.25 | wc -l
'''
