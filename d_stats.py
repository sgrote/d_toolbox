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
	parser = argparse.ArgumentParser(description='Takes pop1-pop4 and pop-info from files, vcf from stdin and calculates blockwise ABBA and BABA for all combinations of pop1-pop4. Writes to file "out_blocks" containing lines with the 4 pops and ABBA, BABA for every block.', usage='zcat file.vcf.gz | d_stats.py pop1 pop2 pop3 pop4 -i pop-info -n 20 -t')
	# mandatory
	parser.add_argument("chrom", help="chromosome number, e.g. '21'") 
	parser.add_argument("pop1", help="File with populations like in pop-column of info-file. Used in position one of D-stats. Or (NEW) a comma-separated string of populations - only with option -l, --poplist.")
	parser.add_argument("pop2", help="Like pop1 for position two of D-stats.")
	parser.add_argument("pop3", help="Like pop1 for position three of D-stats.")
	parser.add_argument("pop4", help="Like pop1 for position four of D-stats.")
	parser.add_argument("-i", "--info", required=True, help="mandatory: csv-file with population metadata for pop1-pop4 (might contain more pops). Needs columns [sample; population; sex]; with sample being identical to the name in the vcf.")
	# blocksize
	parser.add_argument("-k", "--blocksize", type=int, default=200, help="Size of blocks in kb.")

	# optional
	parser.add_argument("-t", "--transver", action="store_true", help="Transversions only.")
	parser.add_argument("-c" ,"--centro", help="optional: bed-file containing chr,start,stop for the centromere of the chromosome. If defined the centromere will be excluded.")
	parser.add_argument("-s", "--sitesfile", choices=["sites", "full"], help="optional: if 'sites' write a 'sites'-file with [block, pos] for every informative position. If 'full' the 'sites'-file also contains all derived allele freqs at informative positions for the populations provided.")

	args = parser.parse_args()

	# vcf from stdin
	if sys.stdin.isatty():
		sys.exit("Error: Needs vcf from stdin, usage: 'zcat file.vcf.gz | d_stats.py pop1 pop2 pop3 pop4 -i pop-info -n 20 -t'")
	vcf = sys.stdin

	########

	## get centromere range
	centro_range = D.get_range(args.chrom ,args.centro) if args.centro else None
	#print(centro_range)

	# skip header until sample line
	vcf_header = D.get_sample_header(vcf)
	#print(vcf_header)

	## get populations
	pops = D.get_pops_from_args(args.pop1, args.pop2, args.pop3, args.pop4)

	## get population matches pop1, pop2, pop3, pop4
	pw_pops = D.get_pw_pops(pops)
	#D.print_table(pw_pops)
	
	## get unique pops
	pops = D.get_unique_pops(pw_pops)

	# from info-file and vcf-header: get col-numbers for every relevant pop and sex for every col
	pop_colnums, col_gender = D.get_pop_colnumbers(args.info, vcf_header, pops)
	#print(pop_colnums)
	#print(col_gender)

	sys.exit()

	# compute blockwise ABBA  [ [[abba][baba]] [[abba][baba]] ...] (NEW param-Afr fixed)
	blocks, sites_comp = D.abba_block_sums(vcf, pop_colnums, pw_pops, col_gender, args.blocksize, centro_range, args.transver, args.sitesfile)

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



