#!/usr/bin/python3

'''
Clone of d_stats.py
Calculate AAAA, BAAA and ABAA and sites per 4-pop-combi
separately for transitions and transversions
for error estimation
get pop1-pop4 from file, vcf from stdin
only fixed sites
'''


import sys
import argparse

import d_functions as D



# list of AAAA, BAAA, ABAA counts per pw_pop [pattern, pw_pop]
# [ [AAAA], [BAAA_tv], [ABAA_tv], [BAAA_ti], [ABAA_ti] ]
# CUATION this assumes pure genotype-input
def count_patterns(vcf, pop_colnums, pw_pops):
	pops = set(pw_pops[0]+pw_pops[1]+pw_pops[2]+pw_pops[3])
	n_comp = len(pw_pops[0])
	nbial_county = 0
	sites_count = 0
	# initialize counts [pattern, pw_pop]
	pattern_counts = [[0]*n_comp, [0]*n_comp, [0]*n_comp, [0]*n_comp, [0]*n_comp]
	
	for line in vcf:
		line = line.rstrip().split()
		## bases/genotypes
		ref = line[3]
		alt = line[4]
		## reduce to biallelic 
		if alt not in ["A","C","T","G","."] or ref not in ["A","C","T","G"]:
			nbial_county += 1
			continue
		
		# nr of fully parsed lines
		sites_count += 1

		## get alternative allele freqs for input populations
		p = {}
		for pop in pops:
			p[pop] = D.get_p(line, pop_colnums[pop])
		
		## check that not all None (pop with variant might not be part of D-stats)
		if set(p.values()) == {None}:
			continue
		
		## is it a transition?
		transi = False
		if ("A" in ref+alt and "G" in ref+alt) or ("C" in ref+alt and "T" in ref+alt):
			transi = True

		## is it AAAA, BAAA, ABAA?  
		for i in range(n_comp):
			p1 = p[pw_pops[0][i]]
			p2 = p[pw_pops[1][i]]
			p3 = p[pw_pops[2][i]]
			p4 = p[pw_pops[3][i]]
			# require fixed and complete
			pset = set([p1,p2,p3,p4])
			if not pset.issubset({0,1}):
				continue
			# [ [AAAA], [BAAA_tv], [ABAA_tv], [BAAA_ti], [ABAA_ti] ]
			# AAAA
			if p1 == p2 == p3 == p4:
				pattern_counts[0][i] += 1
			# BAAA
			elif p2 == p3 == p4 != p1:
				#print([p1,p2,p3,p4])
				#print("BAAA")
				if transi:
					pattern_counts[3][i] += 1
				else:
					pattern_counts[1][i] += 1
			# ABAA
			elif p1 == p3 == p4 != p2:
				#print([p1,p2,p3,p4])
				#print("ABAA")
				if transi:
					pattern_counts[4][i] += 1
				else:
					pattern_counts[2][i] += 1
		if sites_count % 200000 == 0:
			print(line[1])
			
	print("Number of skipped non-biallelic sites: %d" % nbial_county)
	return pattern_counts
	



def main():
	parser = argparse.ArgumentParser(description='Takes pop1-pop4 and pop-info from files, vcf from stdin and counts AAAA, BAAA and ABAA patterns for all combinations of pop1-pop4 (transitions and transversions separately). Writes to file "AAAA_BAAA_ABAA_counts" containing lines with the 4 pops and counts.', usage='zcat file.vcf.gz | count_AAAA_BAAA_ABAA.py pop1 pop2 pop3 pop4 -i pop-info -t')
	# mandatory: either single pops:
	parser.add_argument("--pop1", help="File with populations like in pop-column of info-file. Used in position one of D-stats. Or a comma-separated string of populations. Alternatively --pwpops can be defined.")
	parser.add_argument("--pop2", help="Like pop1 for position two of D-stats.")
	parser.add_argument("--pop3", help="Like pop1 for position three of D-stats.")
	parser.add_argument("--pop4", help="Like pop1 for position four of D-stats.")
	# or a csv-file:
	parser.add_argument("--pwpops", help="csv-file with 4 columns named pop1, pop2, pop3, pop4. This is an alternative input to pop1-pop4 and implies --fixedpairs. Additional columns don't matter.")
	parser.add_argument("-i", "--info", help="mandatory for vcf-input: csv-file with population metadata for pop1-pop4 (might contain more pops). Needs columns [sample; population; sex]; with sample being identical to the name in the vcf.")
	parser.add_argument("-f", "--fixedpairs", action="store_true", help="Don't create all possible combinations of pop1-pop4 but just paste them (all need same length then). --pwpops input implies --fixedpairs.")


	args = parser.parse_args()

	# vcf from stdin
	if sys.stdin.isatty():
		sys.exit("Error: Needs vcf from stdin, usage: 'zcat file.vcf.gz | count_AAAA.py pop1 pop2 pop3 pop4 -i pop-info -n 20 -t'")
	vcf = sys.stdin

	########

	## skip header until sample line
	vcf_header = D.get_sample_header(vcf)
	#print(vcf_header)

	## get population-matches
	if args.pop1:
		pops = D.get_pops_from_args(args.pop1, args.pop2, args.pop3, args.pop4)
		if args.fixedpairs and D.check_pw_pops(pops):
			pw_pops = pops  # pop-matches defined by order
		else:
			pw_pops = D.get_pw_pops(pops)  # create all possible combination
	elif args.pwpops:
		pw_pops = D.get_pw_pops_from_file(args.pwpops)  # pop-matches defined in file
	
	## get unique pops
	pops = D.get_unique_pops(pw_pops)

	# get col-numbers for every relevant pop (and sex for every col for VCF input)
	pop_colnums, col_gender = D.get_pop_colnumbers(args.info, vcf_header, pops)

	# count AAAA, BAAA, ABAA per pop-combi
	pattern_counts = count_patterns(vcf, pop_colnums, pw_pops)

	# rearrange output and print to file [[pop1][pop2][pop3][block][aaaa][baaa]...]
	counts_out = pw_pops + pattern_counts
	
	with open("AAAA_BAAA_ABAA_counts","w") as out:
		out.write('\t'.join(["pop1","pop2","pop3","pop4","aaaa","baaa_tv","abaa_tv","baaa_ti","abaa_ti"])+"\n")
	D.print_table_to_file(counts_out, "AAAA_BAAA_ABAA_counts", mode="a")





if __name__ == "__main__":
	main()



