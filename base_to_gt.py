#!/usr/bin/python3

'''
Take a vcf-file from stdin that has bases in LAST columns instead of genotype, convert bases to fake genotype, modify ALT if needed. Write to stdout.
(use case: Kays BamSNPAddMaf was used to add ape-bases to vcf).
CAUTION: Base-columns are assumed to have no sample-name in header. It also removes extra char that is added to header by Kays script (visible with less-command).
Usage: zcat vcf_with_bases.vcf.gz | base_to_gt.py "sampleA sampleB"
'''

import sys
import argparse

# TODO: check if header has strange characters - if file was not produced by bamSNPAddMaf it will have last 2 chars deleted from header
# maybe TODO: allow bases to be in custom columns, not only in last ones

### helper
def check_alts(alt, base):
	''' check if one of the ALT alleles matches base '''
	alt_alls = alt.split(",")
	for i, b in enumerate(alt_alls):
		if b == base: 
			gt = "/".join(2*str(i+1))
			return alt, gt
	alt = alt + "," + base
	gt = "/".join(2*str(len(alt_alls)+1))
	return alt, gt

''' test
check_alts("A","A") == ('A', '1/1')
check_alts("A","G") == ('A,G', '2/2')
check_alts("A,G","A") == ('A,G', '1/1')
check_alts("A,G","G") == ('A,G', '2/2')
check_alts("A,G","C") == ('A,G,C', '3/3')
check_alts("A,G,C","T") == ('A,G,C,T', '4/4')

'''


def convert_base(ref, alt, base):
	''' convert one base to GT and update ALT '''
	base = base.upper()
	if base == "N" or base == ".":
		gt = "./."			
	elif base == ref:
		gt = "0/0"
	elif alt == ".":
		gt = "1/1"
		alt = base
	else:
		alt, gt = check_alts(alt, base)
	return alt, gt

''' test
convert_base("A", "C", "A") == ('C', '0/0')
convert_base("A", "C", "C") == ('C', '1/1')
convert_base("A", ".", "G") == ('G', '1/1')
convert_base("A", "C", "G") == ('C,G', '2/2')
convert_base("A", "C", "N") == ('C', './.')
convert_base("A", "C", ".") == ('C', './.')
convert_base("A", "C,G", "T") == ('C,G,T', '3/3')

'''


def convert_bases(ref, alt, bases):
	''' convert list of bases to list of GT and update ALT for every base '''
	gts = []
	for base in bases:
		alt, gt = convert_base(ref, alt, base)
		gts.append(gt)
	return alt, gts

''' test
convert_bases("A", "C", ["A","A"]) == ('C', ['0/0', '0/0'])
convert_bases("A", "C", ["C","N","A"]) == ('C', ['1/1', './.', '0/0'])
convert_bases("A", ".", [".","G"]) == ('G', ['./.', '1/1'])
convert_bases("A", "C", ["G","A","T","C","G"]) == ('C,G,T', ['2/2', '0/0', '3/3', '1/1', '2/2'])
convert_bases("A", "C,G", ["T","C",".","A","T"]) == ('C,G,T', ['3/3', '1/1', './.', '0/0', '3/3'])

'''

###########################################################


def main():
	parser = argparse.ArgumentParser(description='Take a vcf-file from stdin that has bases in LAST columns instead of genotype, convert bases to fake genotype, modify ALT if needed. Write to stdout. CAUTION:  use case: Kays BamSNPAddMaf was used to add ape-bases to vcf, so base-columns are assumed to have no sample-name in header. It also deletes last 2 chars in header lines to remove strange sign that is added by Kays script (visible with less-command).', usage='zcat vcf_with_bases.vcf.gz | base_to_gt.py "sampleA sampleB"')
	parser.add_argument("new_names", help="names of samples with bases instead of genotypes, all in one string, ws-separated")

	args = parser.parse_args()

	if sys.stdin.isatty():
		sys.exit("Error: Needs vcf from stdin, e.g. 'zcat vcf_with_bases.vcf.gz | base_to_gt.py | 'sampleA sampleB'")
	instream = sys.stdin

	# how many columns with bases? convert new names to tab-separated string
	args.new_names = args.new_names.split()
	n = len(args.new_names)
	args.new_names = "\t".join(args.new_names)

	# header
	vcf = instream.readline()
	while vcf[:6] != "#CHROM":
		sys.stdout.write(vcf[:(-3)*n]+"\n") # remove strange last character that addMAf adds to header
		vcf = instream.readline()
	last_header = vcf[:(-2)*n] + args.new_names + "\n"
	sys.stdout.write(last_header)

	# body
	for line in instream:
		line = line.strip().split()
		bases = line[-n:]
		ref = line[3]
		alt = line[4]
		# go through new bases, adjust ALT and convert to genotype
		new_alt, new_gts = convert_bases(ref, alt, bases)
		out = line[:4] + [new_alt] + line[5:-n] + new_gts
		sys.stdout.write("\t".join(out) + "\n")


if __name__ == "__main__":
    main()

	

	
