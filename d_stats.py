#!/usr/bin/python3

'''
Calculate blockwise ABBA and BABA and write to file "out_blocks"
containing lines with the 4 pops and ABBA, BABA for every block
get pop1-pop4 from file, vcf from stdin
'''


import sys
import argparse

import d_functions as D


def main():
	parser = argparse.ArgumentParser(description='Takes pop1-pop4 and pop-info from files, vcf from stdin and calculates blockwise ABBA and BABA for all combinations of pop1-pop4. Writes to file "out_blocks" containing lines with the 4 pops and ABBA, BABA for every block.', usage='zcat file.vcf.gz | d_stats.py pop1 pop2 pop3 pop4 -i pop-info -t')
	# mandatory: either single pops:
	parser.add_argument("--pop1", help="File with populations like in pop-column of info-file. Used in position one of D-stats. Or a comma-separated string of populations. Alternatively --pwpops can be defined.")
	parser.add_argument("--pop2", help="Like pop1 for position two of D-stats.")
	parser.add_argument("--pop3", help="Like pop1 for position three of D-stats.")
	parser.add_argument("--pop4", help="Like pop1 for position four of D-stats.")
	# or a csv-file:
	parser.add_argument("--pwpops", help="csv-file with 4 columns named pop1, pop2, pop3, pop4. This is an alternative input to pop1-pop4 and implies --fixedpairs. Additional columns don't matter.")

	# optional
	parser.add_argument("-i", "--info", help="mandatory for vcf-input: csv-file with population metadata for pop1-pop4 (might contain more pops). Needs columns [sample; population; sex]; with sample being identical to the name in the vcf.")
	parser.add_argument("-k", "--blocksize", type=int, default=5000, help="Size of blocks in kb.")
	parser.add_argument("-t", "--transver", action="store_true", help="Transversions only.")
	parser.add_argument("-s", "--sitesfile", choices=["sites", "full"], help="optional: if 'sites' write a 'sites'-file with [block, pos] for every informative position. If 'full' the 'sites'-file also contains all derived allele freqs at informative positions for the populations provided.")
	parser.add_argument("-f", "--fixedpairs", action="store_true", help="Don't create all possible combinations of pop1-pop4 but just paste them (all need same length then). --pwpops input implies --fixedpairs.")
	parser.add_argument("-a", "--afinput", action="store_true", help="Input is already vcf-like file with allele-frequencies per population.")
	parser.add_argument("--centro", help="optional: bed-file containing chr,start,stop for the centromere of the chromosome, which will be excluded. Needs --chrom defined.")
	parser.add_argument("--chrom", help="chromosome number, e.g. '21', only needed if --centro is defined")

	args = parser.parse_args()

	# vcf from stdin
	if sys.stdin.isatty():
		sys.exit("Error: Needs vcf from stdin, usage: 'zcat file.vcf.gz | d_stats.py pop1 pop2 pop3 pop4 -i pop-info -n 20 -t'")
	vcf = sys.stdin

	########

	## get centromere range
	centro_range = D.get_range(args.chrom, args.centro) if args.centro else None
	#print(centro_range)

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
	#print("pop-matches:")
	#D.print_table(pw_pops)
	
	## get unique pops
	pops = D.get_unique_pops(pw_pops)

	# get col-numbers for every relevant pop (and sex for every col for VCF input)
	if args.afinput:
		pop_colnums = D.get_sample_colnumber(vcf_header, pops)
		col_gender = None
	else:
		pop_colnums, col_gender = D.get_pop_colnumbers(args.info, vcf_header, pops)

	# compute blockwise ABBA  [ [[abba][baba]] [[abba][baba]] ...]
	blocks, sites_comp = D.abba_block_sums(vcf, pop_colnums, pw_pops, col_gender, args.blocksize, centro_range, args.transver, args.sitesfile, args.afinput)

	# rearrange output and print to file [[pop1][pop2][pop3][block][abba][baba]]
	blocks_out = D.rearrange_blocks(pw_pops, blocks)
	with open("out_blocks","w") as out:
		out.write('\t'.join(["pop1","pop2","pop3","pop4","block","abba","baba"])+"\n")
	D.print_table_to_file(blocks_out, "out_blocks", mode="a")

	# print sites per comparison to file
	pw_pops.append(sites_comp)
	with open("sites_comp","w") as out:
		out.write('\t'.join(["pop1","pop2","pop3","pop4","n_sites"])+"\n")
	D.print_table_to_file(pw_pops, "sites_comp", mode="a")




if __name__ == "__main__":
	main()



