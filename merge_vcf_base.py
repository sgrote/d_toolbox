#!/usr/bin/python3

''' 
Add fake genotype from [chr,pos,base] file to vcf
only for sites present in vcf (ok for use case 'add random-sampled base to archaics+invariable')
(TODO: modify random sampling script to include REF base)

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
	parser.add_argument("--fill2", action="store_true", help="Fill base_file with REF 0/0 for missing positions, else ./. (this will rarely be useful, only when base_file contains only non-REF bases)")
	parser.add_argument("--require2", action="store_true", help="Restrict to positions present in base_file")
	parser.add_argument("--keep_miss", action="store_true", help="Keep vcf lines that were filtered out as missing data ./. even if base is not present. (If filter fails but base present line is always kept.)")
	parser.add_argument("--filter1", help="Expression to filter vcf, e.g. 'QUAL > 0 and AC > 10'. Acts on QUAL, FILTER and INFO fields. CAUTION: needs whitespaces around keywords.")
	parser.add_argument("--filter_ind1", help="Expression to filter individual genotypes of vcf, e.g. 'GF >= 0 or GT == 1/1'. Acts on keywords in FORMAT. Not passing genotypes will be replaced by './.' . CAUTION: needs whitespaces around keywords.")
	
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

	# new genotype if pos is not present in base_file
	if args.fill2:
		fill_base_gt = "0/0"
	else:
		fill_base_gt = "./."

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
				# -> if base not required and vcf passes filter: add base 0/0 or ./. GT and read next vcf line
				if (len(base) == 0) or (int(vcf[1]) < int(base[1])):
					if not args.require2:
						filtered_vcf = F.combi_filter(vcf, args.filter1, args.filter_ind1, args.gt_only, args.var, args.keep_miss)
						if filtered_vcf:
							out = filtered_vcf + [fill_base_gt]
						elif args.keep_miss:
							out = vcf[:9] + ["./."] * (nvcf) + [fill_base_gt]
					vcf = vcf_in.readline().split()
				## pos not vcf --> skip: since REF at this pos is not known
				elif int(base[1]) < int(vcf[1]):
					base = bases.readline().split()
				## positions are present in both files
				elif vcf[1] == base[1]:
					new_alt, base_gt = B.convert_base(vcf[3], vcf[4], base[2])
					if not (args.var and new_alt == "."):
						# apply line-filter and genotype filter (replace non-passing genotypes with ./.)
						# TODO: if ALT is filtered out, line is still printed despite var=True (egal...)
						# (but can't GT-check merged line since FORMAT might differ)
						filtered_vcf = F.combi_filter(vcf, args.filter1, args.filter_ind1, args.gt_only, False, True) # var=False since base might be only ALT; keep_miss=T -> never None
						out = filtered_vcf[:4] + [new_alt] + filtered_vcf[5:] + [base_gt]
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
VCF=/mnt/sequencedb/gendivdata/2_genotypes/mergedArchModernApes/merged_ancient_apes_sgdp1_chr21.vcf.gz
BASE=/mnt/scratch/steffi/D/random_bases/Forbes_Quarry/FQDeam/ForbesDeam_chr21.tab.gz
MERGE=/mnt/expressions/steffi/D/d_toolbox/merge_vcf_base.py
zcat $VCF | $MERGE $BASE ForbesDeam | less -S
zcat $VCF | $MERGE $BASE ForbesDeam --require2 | less -S
zcat $VCF | $MERGE $BASE ForbesDeam --require2 --var | less -S
zcat $VCF | $MERGE $BASE ForbesDeam --filter1 "QUAL > 100" --var | less -S
zcat $VCF | $MERGE $BASE ForbesDeam --filter1 "QUAL > 100" --keep_miss | less -S

'''

	
