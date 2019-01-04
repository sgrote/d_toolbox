#!/usr/bin/python3

''' 
Filter one vcf file with header, from STDIN

Remove lines that are only variable due to alleles in a group of samples
Mask genotypes that are private to the group in lines that are kept due to other variation  

Writes to STDOUT
'''

import argparse
import sys

import d_functions as D
import merge_vcf_GTs as M


# TODO: maybe add option to allow lines where allele is present at least x times inside the group

def main():
	parser = argparse.ArgumentParser(description='Filter one vcf file, from STDIN removing lines that are only variable due to alleles in a group of samples. Mask genotypes that are private to the group in lines that are kept due to other variation Writes to STDOUT.', usage='zcat file1.vcf.gz | remove_group_singletons.py group_file')
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
	group_index = D.get_sample_colnumber(sample_header, samples).values()
	
	# all other colnums with genotypes
	outside_index = [i for i in range(9, len(sample_header)) if i not in group_index]

	county_invar = 0
	county_private = 0
	county_masked = 0

	for line in vcf:
		line = line.rstrip().split()
		try:
			# no Alt
			if line[4] == ".":
				if args.var:
					county_invar += 1
				else:
					sys.stdout.write("\t".join(line) + "\n")
			else:
				# invariable outside the group and --var? -> remove
				outside_alleles = D.unique_variants(line, outside_index) # e.g [0], [0,2]
				if args.var and len(outside_alleles) < 2:
					# ALT allele should come from group (despite cases of all-ALT, i.e. county_invar - egal)
					county_private += 1 
					continue
				
				# variable outside group or not --var?  --> check for singletons and mask
				inside_alleles = D.unique_variants(line, group_index)
				group_single = [i for i in inside_alleles if i not in outside_alleles]
				if len(group_single) > 0:
					# mask
					for s in group_index:
						for g in group_single:
							line[s] = line[s].replace(g,".")
					# update ALT
					line = M.check_complete_alt(line)
					county_masked += 1
				sys.stdout.write("\t".join(line) + "\n")
				
		except (IndexError, ValueError) as errore:
			sys.stderr.write(str(errore) + "\n")
			sys.stderr.write("\t".join(line[:12]) + "\n")
			sys.exit()
			
	sys.stderr.write("Nr of removed invariable sites: " + str(county_invar) + "\n")
	sys.stderr.write("Nr of removed group singleton sites: " + str(county_private) + "\n")
	sys.stderr.write("Nr of masked group singleton sites: " + str(county_masked) + "\n")


if __name__ == "__main__":
	main()
