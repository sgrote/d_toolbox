#!/usr/bin/python3

'''
Calculate allele freqs per population defined in pops and info file
optionally restrict to transversions
similar to vcf-file to also use it with d_stats.py
header: #CHROM  POS  FLAG  REF  ALT  AltaiNeandertal  Mbuti ...
write to stdout
take vcf from stdin
'''

import sys
import argparse
import gzip

import d_functions as D


def main():
	parser = argparse.ArgumentParser(description='Takes pops and pop-info from files, vcf from stdin and calculates allele freqs per population. Writes to stdout [#CHROM, POS, FLAG, REF, ALT, freq_pop1, freq_pop2, ...]. Only for biallelic sites, omits sites where all freqs are in {0, None}.', usage='zcat chr21.vcf.gz | allele_freqs.py pops -i pop-info | bgzip > alle_freqs_chr21.tab.gz')
	# mandatory
	parser.add_argument("pops", help="File with populations like in pop-column of info-file.")
	parser.add_argument("-i", "--info", required=True, help="mandatory: csv-file with population metadata for pops (might contain more pops). Needs columns [sample; population; sex]; with sample being identical to the name in the vcf.")
	# optional
	parser.add_argument("-t", "--transver", action="store_true", help="Transversions only.")
	parser.add_argument("-d", "--digits", type=int, help="Round allele frequencies to -d digits.")
	parser.add_argument("-c", "--counts", action="store_true", help="Output 'ALT/total' allele counts instead of frequencies, e.g '3/6' instead of '0.5'")
	parser.add_argument("-s", "--store_failed", help="Optional file to write lines to that were filtered out, e.g. because they are not biallelic. For sites with more than 1 ALT allele, this will also have multiple allele frequencies per population.")

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
	
	# open file for filtered-out sites
	fail_handle = open_if(args.store_failed)
	if fail_handle:
		fail_handle.write('\t'.join(["#CHROM","POS","FLAG","REF","ALT"] + list(pops)) + "\n")
	
	# compute allele freqs and write to stdout
	allele_freqs(vcf, pop_colnums, pops, args.transver, args.digits, args.counts, fail_handle)
	
	# close file for filtered-out sites
	if fail_handle:
		fail_handle.close()


def open_if(outfile):
	if not outfile:
		return None
	if not outfile.endswith(".gz"):
		outfile = outfile + ".gz"
	out = gzip.open(outfile, "wt")
	return out


def add_flag(old_flag, new_flag):
	if old_flag == ".":
		return new_flag
	if new_flag in old_flag:
		return old_flag
	else:
		return old_flag + "," + new_flag


def flag(fields):
	'''
	evaluate fields and return flag
	'''
	flag = "."
	# check REF
	ref = fields[3]
	if len(ref) > 1 or ref=="-":
		flag = add_flag(flag, "indel")
	elif ref not in ["A","C","T","G","-"]:
		flag = add_flag(flag, "refbase")
	# check ALTs
	if fields[4] == ".":
		flag = add_flag(flag, "invar")
		return flag
	alts = fields[4].split(",")
	if len(alts) > 1:
		flag = add_flag(flag, "multial")
	if any([len(a) > 1 for a in alts]) or "-" in alts:
		flag = add_flag(flag, "indel")
	onebases = [a for a in alts if len(a) == 1]
	if any([o not in ["A","C","T","G","-"] for o in onebases]):
		flag = add_flag(flag, "altbase")
	# any transition among single base REF-ALT?
	if len(ref) == 1:
		bases = onebases + [ref]
		if ("A" in bases and "G" in bases) or ("C" in bases and "T" in bases):
			flag = add_flag(flag, "transi")
	return flag

''' test
flag(["21", "1234", ".", "A", "T"]) == "."
flag(["21", "1234", ".", "A", "."]) == "invar"
flag(["21", "1234", ".", "N", "."]) == "refbase,invar"
flag(["21", "1234", ".", "TTT", "."]) == "indel,invar"
flag(["21", "1234", ".", "TTT", "C"]) == "indel"
flag(["21", "1234", ".", "TTT", "CC"]) == "indel"
flag(["21", "1234", ".", "A", "G"]) == "transi"
flag(["21", "1234", ".", "T", "C"]) == "transi"
flag(["21", "1234", ".", "T", "CC"]) == "indel"
flag(["21", "1234", ".", "T", "C,TTT"]) == "multial,indel,transi"
flag(["21", "1234", ".", "T", "C,TTT,N"]) == "multial,indel,altbase,transi"

'''


def get_p_multi(vcf_line, pop, n_alts, digits=None, counts=False):
	'''
	clone of D.get_p that also works with multiallelic sites
	alternative allele-freqs for one population
	in: vcf = vcf-line as list
		pop = [colnumbers] for the population
		n_alts = nr of alternative alleles
		digits = digits to round to
		counts = return ("ALT1,ALT2/total") counts as string instead of the fraction
	out: freq1,freq2 as string for all ALTs
	'''
	
	## pre-check for all-REF or all missing
	genotypes = [vcf_line[colnumber] for colnumber in pop]
	if not counts and all([g in ["0/0", "0|0", "0"] for g in genotypes]):
		return "0.0"
	if all([g in ["./.", ".|.", "."] for g in genotypes]):
		return None
	
	# count alternative alleles and total number of alleles for one population
	ac = [0] * n_alts
	sum_alleles = 0
	for i in range(len(pop)):
		gt = genotypes[i].replace('/', '|')
		alleles = gt.split("|")
		for a in alleles:
			if a == ".":
				continue
			if a == "0":
				sum_alleles += 1
			else:
				# "1" = ac[0] etc
				a = int(a)
				ac[a-1] += 1
				sum_alleles += 1

	# compute frequency
	if counts:
		daf = ",".join(map(str, ac)) + "/" + str(sum_alleles)
	else:
		if digits:
			daf = [round(a*1.0 / sum_alleles, digits) for a in ac]
		else:
			daf = [a*1.0 / sum_alleles for a in ac]
		daf = ",".join(map(str, daf))
	return daf


''' test
vcf_line = ['21','148','.','C','T','0','.','.','GT','0/1','1/1','./0','./.','1|.','1|1','.|.','0/0','0|0']
get_p_multi(vcf_line, [11], 1) == '0.0'
get_p_multi(vcf_line, [9,10], 1) == '0.75'
get_p_multi(vcf_line, [9,10,11], 1) == '0.6'
vcf_line = ['21','148','.','C','T,GT','0','.','.','GT','0/1','1/2','./0','./.','1|.','2|2']
get_p_multi(vcf_line, [9,10,11], 2) == '0.4,0.2'
get_p_multi(vcf_line, [9,10,11], 2, counts=True) == '2,1/5'
get_p_multi(vcf_line, [13,14], 2, digits=3) == '0.333,0.667'
get_p_multi(vcf_line, [12], 2) == None
get_p_multi(vcf_line, [11,12], 2) == '0.0,0.0'
vcf_line = ['21','148','.','C','T,GT,A','0','.','.','GT','3/1','1/2','./0','./.','1|.','2|2']
get_p_multi(vcf_line, [9,10], 3) == '0.5,0.25,0.25'

'''


# CAUTION this assumes pure genotype-input
def allele_freqs(vcf, pop_colnums, pops, transver=False, digits=None, counts=False, fail_handle=None):
	sys.stdout.write('\t'.join(["#CHROM","POS","FLAG","REF","ALT"] + list(pops)) + "\n")
	line_count = 0
	trans_county = 0
	nbial_county = 0
	# unique colnumbers of input and non-input samples
	pop_cols, non_pop_cols = D.unique_pop_colnums(pop_colnums, len(pops))
	for line in vcf:
		try:
			fields = line.rstrip().split()
			# modify fields if more than 2 states in ALT but only 2 in samples
			# returns None or an updated biallelic fields
			fields = D.check_bial(fields, pop_cols, non_pop_cols)
			if not fields:
				# if no pass get unmodified fields again, will be flagged as multi-allelic
				if fail_handle:
					fields = line.rstrip().split()
				else:
					continue
			# get flag
			fields[2] = flag(fields)
			if fields[2] != "." and not fail_handle:
				continue
			# get freqs
			n_alts = len(fields[4].split(","))
			freqs = []
			for pop in pops:
				freqs.append(get_p_multi(fields, pop_colnums[pop], n_alts, digits, counts))
			# replace 'None' with '.'
			freqs = ["." if p == None else p for p in freqs]
			# write to stdout or file with fails
			if fields[2] == ".":
				sys.stdout.write("\t".join(fields[:5]) + "\t" + "\t".join(map(str,freqs)) + "\n")
			else:
				fail_handle.write("\t".join(fields[:5]) + "\t" + "\t".join(map(str,freqs)) + "\n")
			
		except(IndexError, ValueError) as errore:
			sys.stderr.write(str(errore) + "\n")
			sys.stderr.write("\t".join(fields) + "\n")
			sys.exit()





if __name__ == "__main__":
	main()





