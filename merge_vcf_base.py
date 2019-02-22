#!/usr/bin/python3

''' 
Add fake genotype from [chr,pos,base] file to vcf
only for sites present in vcf (ok for use case 'add random-sampled base to archaics+invariable')

'''
	
import sys
import argparse
import gzip

import vcf_line_filter as F
import base_to_gt as B

def main():
	parser = argparse.ArgumentParser(description='Merge a vcf file from STDIN with a file containing [chr | pos | base] (compressed or uncompressed). Optional filtering of vcf. Header is taken from vcf; positions not in vcf are skipped. Writes to STDOUT.', usage='zcat file1.vcf.gz | merge_vcf_base.py base_file.gz newSample --filter1 "QUAL > 3" --require2')
	parser.add_argument("base_file", help="File with chrom, pos, base; 1-based")
	parser.add_argument("base_name", help="Name of new column in the merged vcf-header")
	parser.add_argument("--var", action="store_true", help="Restrict to variable sites")
	parser.add_argument("--gt_only", action="store_true", help="Restrict to genotypes in sample columns")
	parser.add_argument("--require2", action="store_true", help="Restrict to positions present in base_file")
	
	args = parser.parse_args()

	# check instream
	if sys.stdin.isatty():
		sys.exit("Error: Needs vcf from stdin, e.g. 'zcat file1.vcf.gz | merge_vcf_base.py base_file.gz newSample'")
	vcf_in = sys.stdin
	
	# print header, add new sample to line with donor-ids
	vcf = vcf_in.readline()
	while vcf[:6] != "#CHROM":
		sys.stdout.write(vcf)
		vcf = vcf_in.readline()
	sys.stdout.write("\t".join([vcf.rstrip()] + [args.base_name]) + "\n")

	# open base-file
	if args.base_file.endswith('.gz'):
	    opener = gzip.open
	else:
	    opener = open
	with opener(args.base_file, 'rt') as bases:
		vcf = vcf_in.readline().split()
		base = bases.readline().split()
		if vcf[0] != base[0]:
			sys.exit("Error: Chromosomes do not match: vcf-chrom={}, base-chrom={}.".format(vcf[0], base[0]))
		## number of donors
		n_vcf = len(vcf)-9
		### iterate over lines until vcf ends (pos not in vcf are ignored)
		while not(len(vcf) == 0):
			## for every line: check if line passes and filter on individual genotypes
			try:
				out = None
				# pos not in base, including the case that base file has reached end
				# -> if base not required: add base ./. GT and read next vcf line
				if (len(base) == 0) or (int(vcf[1]) < int(base[1])):
					if not args.require2:
						out = vcf + ["./."]
					vcf = vcf_in.readline().split()
				## pos not vcf --> skip: since REF at this pos is not known
				elif int(base[1]) < int(vcf[1]):
					base = bases.readline().split()
				## positions are present in both files
				elif vcf[1] == base[1]:
					new_alt, base_gt = B.convert_base(vcf[3], vcf[4], base[2])
					if not (args.var and new_alt == "."):
						out = vcf[:4] + [new_alt] + vcf[5:] + [base_gt]
					# read next lines from both files (also if not (vcf_pass or base_pass))
					vcf = vcf_in.readline().split()
					base = bases.readline().split()
				# print merged line if any
				if out:
					sys.stdout.write("\t".join(out) + "\n")
			except (IndexError, ValueError) as errore:
				sys.stderr.write(str(errore) + "\n")
				sys.stderr.write("\t".join(vcf) + "\n")
				sys.stderr.write("\t".join (base) + "\n")
				sys.exit()


if __name__ == "__main__":
	main()

''' test
VCF=/mnt/sequencedb/gendivdata/2_genotypes/giantVcfs/merged_all_sites_arch_apes_sgdp1_g1000_chr21.vcf.gz
BASE=/mnt/scratch/steffi/D/random_bases/Forbes_Quarry/forbes_deam_chr21.tab.gz
MERGE=/mnt/expressions/steffi/D/d_toolbox/merge_vcf_base.py
zcat $VCF | cut -f-15 | $MERGE $BASE ForbesDeam | less -S
zcat $VCF | cut -f-15 | $MERGE $BASE ForbesDeam --require2 | less -S
zcat $VCF | cut -f-15 | $MERGE $BASE ForbesDeam --require2 --var | less -S

'''

	
