#!/usr/bin/python3

'''
Calculate pairwise divergences between populations defined in pops and info file
optionally restrict to transversions
write to file "out_diverge" [pop1, pop2, sum_match, sum_mismatch]
take vcf from stdin
'''

import sys
import argparse
import re

import d_functions as D
import allele_freqs as A


def main():
	parser = argparse.ArgumentParser(description='Calculate pariwise divergences between populations defined in pops and info file. Writes to file "out_diverge" [pop1, pop2, sum_match, sum_mismatch].', usage='zcat chr21.vcf.gz | pw_diverge.py pops -i pop-info')
	# mandatory
	parser.add_argument("pops", help="File with populations like in pop-column of info-file.")
	parser.add_argument("-i", "--info", required=True, help="mandatory: csv-file with population metadata for pops (might contain more pops). Needs columns [sample; population; sex]; with sample being identical to the name in the vcf.")
	# optional
	parser.add_argument("--transver", action="store_true", help="Transversions only.")
	parser.add_argument("--homo", action="store_true", help="All input genotypes are homozygous (like for fake genotypes from randomly sampled alleles) - faster computation.")

	args = parser.parse_args()
	print(args)
	
	# vcf from stdin
	if sys.stdin.isatty():
		sys.exit("Error: Needs vcf from stdin.")
	vcf = sys.stdin

	########

	# skip header until sample line
	vcf_header = D.get_sample_header(vcf)

	## get list of populations
	print("Read pops file...")
	pops = A.get_pops_from_file(args.pops)
	
	# from info-file and vcf-header:
	# get col-numbers for every relevant pop and sex for every col (though sex not needed here)
	print("Read info file...")
	pop_colnums, col_gender = D.get_pop_colnumbers(args.info, vcf_header, pops)
	#print(pop_colnums)
	#print(col_gender)

	print("Compute pairwise divergences...")
	out_diverge = pw_diverge(vcf, pop_colnums, pops, args.transver, args.homo)
	
	with open("out_diverge","w") as out:
		out.write('\t'.join(["pop1","pop2","sum_match","sum_mismatch"]) + "\n")
	D.print_table_to_file(out_diverge, "out_diverge", mode="a")


def get_pw_pops(pop_list):
	'''
	in: list [pop1, pop2, ..., popn]
	out: list of 2 lists [[pop1], [pop2]] for all pairwise matches
	'''
	pw_pops = [[],[]]
	for i in range(len(pop_list)):
		for j in range(i+1, len(pop_list)):
			pw_pops[0].append(pop_list[i])
			pw_pops[1].append(pop_list[j])
	return pw_pops

''' test
pops=["A","B","C","D"]
get_pw_pops(pops)
'''


def get_pop_genotypes(pop_colnums, vcf_line):
	'''
	in: dict{pop:[colnumbers]}
		vcf_line [chr, pos, ..., gt_sample1, gt_sample2]
	out: dict{pop:[genotypes]}
	'''
	pop_genotypes = {}
	for pop in pop_colnums.keys():
		pop_genotypes[pop] = [vcf_line[i] for i in pop_colnums[pop]]
	return pop_genotypes

''' test
pop_colnums = {"A":[5,6], "B":[8,4,7]}
vcf_line = ["21", "23536", "C", "T", "./.", "0/1", "0/0", "1/0", "./1"]
get_pop_genotypes(pop_colnums, vcf_line)
'''


def gts_to_alleles(gts):
	'''
	in: list of genoytpes
	out: list of alleles
	'''
	gts_string = "|".join(gts)
	alleles = re.split("[/|]", gts_string)
	return alleles

''' test
gts_to_alleles(['0/1', '1|1', '.', '1/.']) == ['0', '1', '1', '1', '.', '1', '.']
'''


def comp_one_homo (gt1, gt2):
	'''
	compare divergence of one pop-match with all homozygous genotypes
	'''
	match = 0
	mismatch = 0
	# remove "./." and ".|." before comparing
	gt1_called = [gt for gt in gt1 if gt not in ["./.", ".|."]]
	gt2_called = [gt for gt in gt2 if gt not in ["./.", ".|."]]
	for g1 in gt1_called:
		g1 = g1.replace("|","/")
		for g2 in gt2_called:
			g2 = g2.replace("|","/")
			# check for match
			if g1 == g2:
				match += 1
			else:
				mismatch += 1
	return match, mismatch

''' test
comp_one_homo(['1/1', '0|0'], ['./.', '0|0']) == (1, 1)
comp_one_homo(['1/1'], ['./.', '0|0']) == (0, 1)
comp_one_homo(['0|0', '0/0', '.|.'], ['./.', '0|0']) == (2, 0)
'''

def comp_one(gt1, gt2):
	'''
	compare divergence of one pop-match between single alleles
	'''
	match = 0
	mismatch = 0
	alleles1 = gts_to_alleles(gt1)
	alleles2 = gts_to_alleles(gt2)
	# remove "." before comparing
	a1_called = [a for a in alleles1 if a != "."]
	a2_called = [a for a in alleles2 if a != "."]
	for a1 in a1_called:
		for a2 in a2_called:
			# check for match
			if a1 == a2:
				match += 1
			else:
				mismatch += 1
	
	
	return match, mismatch

''' test
comp_one(['1/1', '0|0'], ['./.', '0|0']) == (4, 4)
comp_one(['1/1'], ['./.', '0|0']) == (0, 4)
comp_one(['0|0', '0/0', '.|.'], ['./.', '0|0']) == (8, 0)
comp_one(['0|1', '0/0', '.|.'], ['./.', '1|0']) == (4, 4)
comp_one(['0|1', '0/0', '.|.'], ['./.', '1|2']) == (1, 7)
comp_one(['0|1'], ['1/1']) == (2, 2)
comp_one(['1|1'], ['./1']) == (2, 0)
'''


def pw_diverge(vcf, pop_colnums, pops, transver=False, homo=False):
	'''
	compute sum of matches and mismatches for each pw_pop
	list of lists [pop1, pop2, sum_match, sum_mismatch]
	'''
	# generate matches
	pw_pops = get_pw_pops(pops)
	D.print_table(pw_pops)
	n_comp = len(pw_pops[0])
	# initialize sums	
	line_count = 0
	trans_county = 0
	nbial_county = 0
	matches = [0] * n_comp
	mismatches = [0] * n_comp
	
	for line in vcf:
		try:
			line = line.rstrip().split()
			## bases/genotypes
			ref = line[3]
			alt = line[4]
			## reduce to biallelic or fixed REF
			# TODO: make biallelic and/or variable optional)
			if alt not in ["A","C","T","G","."] or ref not in ["A","C","T","G"]:
				#print("Skipping non-biallelic:", line[:7])
				nbial_county += 1
				continue
			## reduce to transversions 
			if transver: 
				if ("A" in ref+alt and "G" in ref+alt) or ("C" in ref+alt and "T" in ref+alt):
					#print("Skipping transition:", line[:7])
					trans_county += 1
					continue
			## get list of genotypes for all populations (dict{pop:[genotypes]})
			pop_genotypes = get_pop_genotypes(pop_colnums, line)
			
			## count matches and mismatches
			for i in range(n_comp):
				if homo:
					ma, mi = comp_one_homo(pop_genotypes[pw_pops[0][i]], pop_genotypes[pw_pops[1][i]])
				else:
					ma, mi = comp_one(pop_genotypes[pw_pops[0][i]], pop_genotypes[pw_pops[1][i]])
				# add count for pop-match
				matches[i] += ma
				mismatches[i] += mi
			## line-count to screen
			line_count += 1
			if line_count % 100000 == 0:
				print(line[1])
		except (IndexError, ValueError) as errore:
			print(errore)
			print(line)
			sys.exit()
	print("Number of skipped non-biallelic sites: %d" % nbial_county)
	print("Number of skipped transitions: %d" % trans_county)
	print("Number of valid lines: %d" % line_count)
	## create output
	out = [pw_pops[0], pw_pops[1], matches, mismatches]
	return out



if __name__ == "__main__":
	main()





