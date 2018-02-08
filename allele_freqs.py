#!/usr/bin/python3

'''
Calculate allele freqs per population defined in pops and info file
optionally restrict to transversions
write to file "out_freqs" [chr, pos, ref, alt, freq_pop1, freq_pop2, ...]
take vcf from stdin
'''

import sys
import argparse

import d_functions as D


def main():
	parser = argparse.ArgumentParser(description='Takes pops and pop-info from files, vcf from stdin and calculates allele freqs per population. Writes to file "out_freqs" [chr, pos, ref, alt, freq_pop1, freq_pop2, ...].', usage='zcat chr21.vcf.gz | allele_freqs.py pops -i pop-info')
	# mandatory
	parser.add_argument("pops", help="File with populations like in pop-column of info-file.")
	parser.add_argument("-i", "--info", required=True, help="mandatory: csv-file with population metadata for pops (might contain more pops). Needs columns [sample; population; sex]; with sample being identical to the name in the vcf.")
	# optional
	parser.add_argument("-t", "--transver", action="store_true", help="Transversions only.")

	args = parser.parse_args()

	# vcf from stdin
	if sys.stdin.isatty():
		sys.exit("Error: Needs vcf from stdin, usage: 'zcat chr21.vcf.gz | allele_freqs.py pops -i pop-info'")
	vcf = sys.stdin

	########

	# skip header until sample line
	vcf_header = D.get_sample_header(vcf)

	## get list of populations
	print("Read pops file...")
	pops = get_pops_from_file(args.pops)
	print(pops)
	
	# from info-file and vcf-header: get col-numbers for every relevant pop and sex for every col
	print("Read info file...")
	pop_colnums, col_gender = D.get_pop_colnumbers(args.info, vcf_header, pops)
	for pop in pop_colnums.keys():
		print(pop)
		print(len(pop_colnums[pop]))
	#print(pop_colnums)
	#print(col_gender)

	print("Compute allele freqs...")
	allele_freqs(vcf, pop_colnums, pops, col_gender, args.transver)


def get_pops_from_file(pop_file):
	'''
	in: file with one column containing population names
	out: list [pop1, pop2, ...]
	'''
	pops = []
	with open(pop_file, "r") as p:
		for pop in p:
			pops.append(pop.rstrip())
	return pops

# write "out_freqs" with [chr, pos, ref, alt, freq_pop1, freq_pop2, ...]
# CUATION this assumes pure genotype-input
def allele_freqs(vcf, pop_colnums, pops, col_gender, transver=False):
	with open("out_freqs","w") as out:
		out.write('\t'.join(["chr","pos","ref","alt"] + list(pops)) + "\n")
		line_count = 0
		trans_county = 0
		nbial_county = 0
		for line in vcf:
			try:
				line = line.rstrip().split()
				## bases/genotypes
				ref = line[3]
				alt = line[4]
				## reduce to biallelic 
				if alt not in ["A","C","T","G"] or ref not in ["A","C","T","G"]:
					#print("Skipping non-biallelic:", line[:7])
					nbial_county += 1
					continue
				## reduce to transversions 
				if transver: 
					if ("A" in ref+alt and "G" in ref+alt) or ("C" in ref+alt and "T" in ref+alt):
						#print("Skipping transition:", line[:7])
						trans_county += 1
						continue
				## get alternative allele freqs for input populations
				p = []
				for pop in pops:    
					p.append(D.get_p(line, pop_colnums[pop], col_gender, X=False)) # CAUTION: only for autos
				#print(p)
				## check that not all pop-freqs are None (no genotypes in pop)
				pset = set(p)
				if pset.issubset({0,None}):
					#print("Uninformative site:", pset) 
					continue
				line_count += 1
				if line_count % 50000 == 0:
					print(line_count)
				out.write("\t".join(line[0:2] + line[3:5]) + "\t" + "\t".join(map(str,p)) + "\n")
			except (IndexError, ValueError) as errore:
				print(errore)
				print(line)
	print("Number of skipped non-biallelic sites: %d" % nbial_county)
	print("Number of skipped transitions: %d" % trans_county)
	print("Number of valid lines: %d" % line_count)
	

if __name__ == "__main__":
	main()





