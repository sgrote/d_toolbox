#!/usr/bin/python3

''' 
filter streamed in vcf-like file with bed-file
output to stdout
input stream: (#header), columns: chr, pos, ... (1-based)
input bed-file (compressed or uncompressed): chr, pos, end (0-based, exclusive)

assumes vcf is only one chrom 

'''


import sys
import argparse
import gzip


def main():
	parser = argparse.ArgumentParser(description="filter streamed in vcf-like file with bed-file. Write to stdout. input stream: (#header), columns: chr, pos, ... (1-based). Assumes vcf is only one chrom.", usage="zcat file.vcf.gz | filter_vcf_with_bed.py filter.bed")
	parser.add_argument("bed_file", help="chr, pos, end (0-based, exclusive). Positions to be kept in vcf, or to be removed with option --remove. *.bed or *.bed.gz")
	parser.add_argument("-r", "--remove", action="store_true", help="Positions in bed-files are removed and not kept.")

	args = parser.parse_args()
	
	if sys.stdin.isatty():
		sys.exit("Error: Needs vcf from stdin, e.g. 'zcat file.vcf.gz | filter_vcf_with_bed.py filter.bed'")
	vcf_file = sys.stdin

	## print header
	vcf = vcf_file.readline()
	while vcf[0] == "#":
		sys.stdout.write(vcf)
		vcf = vcf_file.readline()
	vcf = vcf.split()
	
	if args.bed_file.endswith('.gz'):
		opener = gzip.open
	else:
		opener = open
	with opener(args.bed_file, 'rt') as bed_f:

		bed = bed_f.readline().split()

		## check if vcf position is in bed (bed 1-3 = vcf 2-3)
		while not (len(vcf) == 0):

			# avoid BrokenPipeError: [Errno 32] Broken pipe (when bed-file is done but instream not)
			if len(bed) == 0:
				vcf = vcf_file.readline().split()

			# check that chrom matches:
			elif vcf[0] != bed[0]:
				bed = bed_f.readline().split()
			# A) vcf in bed
			elif ((int(vcf[1]) > int(bed[1])) and (int(vcf[1]) <= int(bed[2]))):
				if not args.remove:
					sys.stdout.write("\t".join(vcf) + "\n")
				vcf = vcf_file.readline().split()
			# B) vcf smaller than bed  ->  not in bed
			elif (int(vcf[1]) <= int(bed[1])):
				if args.remove:
					sys.stdout.write("\t".join(vcf) + "\n")
				#print("filtered out:" + "\t".join(vcf[:5]))
				vcf = vcf_file.readline().split()
			# C) vcf larger than bed  ->  skip bed lines until A) or B)
			elif (int(vcf[1]) > int(bed[2])):
				bed = bed_f.readline().split()




if __name__ == "__main__":
	main()
