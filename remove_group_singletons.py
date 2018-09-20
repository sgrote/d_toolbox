#!/usr/bin/python3

''' 
Filter one vcf file with header, from STDIN
Remove lines that are only variable due to alleles in a group of samples
Writes to STDOUT

'''

import argparse
import sys

import d_functions as D

# TODO: maybe add option to allow lines where allele is present at least x times inside the group

def main():
	parser = argparse.ArgumentParser(description='Filter one vcf file, from STDIN removing lines that are only variable due to alleles in a group of samples. Writes to STDOUT.', usage='zcat file1.vcf.gz | remove_group_singletons.py group_file')
	parser.add_argument("group_file", help="text file with one sample name per line, defining the samples in the group.")
	parser.add_argument("--var", action="store_true", help="Restrict to variable sites.")
	
	args = parser.parse_args()

	if sys.stdin.isatty():
		sys.exit('Error: Needs vcf from stdin')
	vcf = sys.stdin


	# get list of group samples
	samples = D.get_pops_from_file(args.group_file)
	
	# vcf header
	line = vcf.readline()
	while line[:6] != "#CHROM":
		sys.stdout.write(line)
		line = vcf.readline()
	# sample-line
	sys.stdout.write(line)
	sample_header = line.rstrip().split()
	
	# get col-numbers for every input sample
	sample_index = D.get_sample_colnumber(sample_header, samples)
	
	# all other colnums with genotypes
	outside_index = [i for i in range(9, len(sample_header)) if i not in sample_index.values()]


	for line in vcf:
		line = line.rstrip().split()
		try:
			# leave invariable
			if line[4] == ".":
				if not args.var:
					sys.stdout.write("\t".join(line) + "\n")
			# for variable sites make sure they are variable outside the group
			else:
				outside_alleles = D.unique_variants(line, outside_index) # unique ALT indices, without "."
				if len(outside_alleles) > 1:
					sys.stdout.write("\t".join(line) + "\n")

		except (IndexError, ValueError) as errore:
			sys.stderr.write(str(errore) + "\n")
			sys.stderr.write("\t".join(line[:12]) + "\n")
			sys.exit()


if __name__ == "__main__":
	main()
