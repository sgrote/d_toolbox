#!/usr/bin/python3

'''
Calculate allele freqs per population defined in pops and info file
optionally restrict to transversions
similar to vcf-file to also use it with d_stats.py
header: #CHROM  POS  ID  REF  ALT  AltaiNeandertal  Mbuti ...
write to stdout
take vcf from stdin
'''

import sys
import argparse

import d_functions as D


def main():
	parser = argparse.ArgumentParser(description='Takes pops and pop-info from files, vcf from stdin and calculates allele freqs per population. Writes to stdout [#CHROM, POS, ID, REF, ALT, freq_pop1, freq_pop2, ...]. Only for biallelic sites, omits sites where all freqs are in {0, None}.', usage='zcat chr21.vcf.gz | allele_freqs.py pops -i pop-info | bgzip > alle_freqs_chr21.tab.gz')
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

	# get list of populations
	pops = get_pops_from_file(args.pops)
	
	# from info-file and vcf-header: get col-numbers for every relevant pop and sex for every col
	pop_colnums, col_gender = D.get_pop_colnumbers(args.info, vcf_header, pops)

	# Compute allele freqs and write to stdout
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

# CAUTION this assumes pure genotype-input
def allele_freqs(vcf, pop_colnums, pops, col_gender, transver=False):
	sys.stdout.write('\t'.join(["#CHROM","POS","ID","REF","ALT"] + list(pops)) + "\n")
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
				continue
			## reduce to transversions 
			if transver: 
				if ("A" in ref+alt and "G" in ref+alt) or ("C" in ref+alt and "T" in ref+alt):
					continue
			## get alternative allele freqs for input populations
			freqs = []
			for pop in pops:    
				freqs.append(D.get_p(line, pop_colnums[pop], col_gender, X=False)) # CAUTION: only for autos
			## check that not all pop-freqs are None (no genotypes in pop)
			pset = set(freqs)
			if pset.issubset({0,None}):
				continue
			# replace 'None' with '.'
			freqs = ["." if p == None else p for p in freqs]
			sys.stdout.write("\t".join(line[:5]) + "\t" + "\t".join(map(str,freqs)) + "\n")
		except(IndexError, ValueError) as errore:
			sys.stderr.write(str(errore) + "\n")
			sys.stderr.write("\t".join(line) + "\n")
			sys.exit()



if __name__ == "__main__":
	main()





