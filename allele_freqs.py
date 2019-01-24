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
	parser.add_argument("-a", "--accu_bial", action="store_true", help="Accurately check for biallelic sites among input pops, other samples can have a third state. (default is checking the ALT-field only).")
	parser.add_argument("-d", "--digits", type=int, help="Round allele frequencies to -d digits.")
	parser.add_argument("-c", "--counts", action="store_true", help="Output 'ALT/total' allele counts instead of frequencies, e.g '3/6' instead of '0.5'")

	args = parser.parse_args()

	# vcf from stdin
	if sys.stdin.isatty():
		sys.exit("Error: Needs vcf from stdin, usage: 'zcat chr21.vcf.gz | allele_freqs.py pops -i pop-info'")
	vcf = sys.stdin
	

	########

	# skip header until sample line
	vcf_header = D.get_sample_header(vcf)

	# get list of populations
	pops = D.get_pops_from_file(args.pops)
	
	# from info-file and vcf-header: get col-numbers for every relevant pop and sex for every col
	pop_colnums, col_gender = D.get_pop_colnumbers(args.info, vcf_header, pops)
	

	# Compute allele freqs and write to stdout
	allele_freqs(vcf, pop_colnums, pops, args.transver, args.accu_bial, args.digits, args.counts)



# CAUTION this assumes pure genotype-input
def allele_freqs(vcf, pop_colnums, pops, transver=False, accu_bial=False, digits=None, counts=False):
	sys.stdout.write('\t'.join(["#CHROM","POS","ID","REF","ALT"] + list(pops)) + "\n")
	line_count = 0
	trans_county = 0
	nbial_county = 0
	# unique colnumbers of input and non-input samples
	pop_cols, non_pop_cols = D.unique_pop_colnums(pop_colnums, len(pops))
	for line in vcf:
		try:
			line = line.rstrip().split()
			## skip invariable
			if line[4] == ".":
				continue
			
			## modify line if more than 2 states in ALT but only 2 in samples
			if accu_bial:
				# returns None or an updated biallelic line
				# to make it work with get_p ALT-allele from samples will be set to 1
				line = D.check_bial(line, pop_cols, non_pop_cols) 
				if not line:
					continue
			
			## check REF and ALT
			## bases/genotypes
			ref = line[3]
			alt = line[4]
			if ref not in ["A","C","T","G"]:
				continue
			if alt not in ["A","C","T","G"]:
				continue
			if transver: 
				if ("A" in ref+alt and "G" in ref+alt) or ("C" in ref+alt and "T" in ref+alt):
					continue
			
			## always check for accu_bial before computing freqs, if it wasnt done yet
			# since alt was already filtered, this will never include GT-updates,
			# only checks if REF and ALT are present in genotypes
			# this is faster then computing freqs first and then remove all 0,None
			if not accu_bial:
				line = D.check_bial(line, pop_cols, non_pop_cols)
				if not line:
					continue
			
			## get alternative allele freqs for input populations
			freqs = []
			for pop in pops:
				freqs.append(D.get_p(line, pop_colnums[pop], digits, counts))
			
			## double-check that not all pop-freqs are None (no genotypes in pop)
			pset = set(freqs)
			if pset.issubset({0,None}) or pset.issubset({1,None}):
				sys.stderr.write("\t".join(line) + "\n")
				sys.stderr.write("How could an invariable site make it through accu_bial check?")
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





