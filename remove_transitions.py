#!/usr/bin/python3

'''
restrict to biallelic transversions given index of REF and ALT column (can be applied to several file-formats)
take file from stdin
write to stdout
'''

import sys
import argparse


def main():
	parser = argparse.ArgumentParser(description='restrict to biallelic transversions given index of REF and ALT column (can be applied to several file-formats). Reads from stdin and writes to stdout.', usage='zcat chr21.vcf.gz | remove_transitions.py > chr21_transversions_only.vcf')
	# mandatory
	parser.add_argument("-r","--ref_index", type=int, default=3, help="Column number of REF-allele, 0-based.")
	parser.add_argument("-a","--alt_index", type=int, default=4, help="Column number of ALT-allele, 0-based.")

	args = parser.parse_args()

	# infile from stdin
	if sys.stdin.isatty():
		sys.exit("Error: Needs file from stdin, usage: 'zcat chr21.vcf.gz | remove_transitions.py > chr21_transversions_only.vcf'")
	infile = sys.stdin

	for line in infile:
		try:
			fields = line.rstrip().split()
			# bases
			ref = fields[args.ref_index]
			alt = fields[args.alt_index]
			# reduce to biallelic 
			if alt not in ["A","C","T","G"]:
				continue
			# reduce to transversions 
			if ("A" in ref+alt and "G" in ref+alt) or ("C" in ref+alt and "T" in ref+alt):
				continue
			# print valid line
			sys.stdout.write(line)
		except (IndexError, ValueError) as errore:
			sys.stderr.write(str(errore) + "\n")
			sys.stderr.write(line + "\n")
			sys.exit()

if __name__ == "__main__":
	main()





