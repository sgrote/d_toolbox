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
	parser.add_argument("pop1", help="File with populations like in pop-column of info-file. Used in position one of D-stats.")
	parser.add_argument("pop2", help="Like pop1 for position two of D-stats.")
	parser.add_argument("pop3", help="Like pop1 for position three of D-stats.")
	parser.add_argument("pop4", help="Like pop1 for position four of D-stats.")
	parser.add_argument("-i", "--info", required=True, help="mandatory: csv-file with population metadata for pop1-pop4 (might contain more pops). Needs columns [sample; population; sex]; with sample being identical to the name in the vcf.")
	# blocksize
	parser.add_argument("-k", "--blocksize", type=int, default=200, help="Size of blocks in kb.")

	# optional
	parser.add_argument("-t", "--transver", action="store_true", help="Transversions only.")
	#parser.add_argument("--startstop", default="/mnt/expressions/steffi/chromSizes/hg19_chrall_ranges.txt", help="optional: bed-file containing chr,start,stop for the boundaries of the chromosome to take into account (ugly but required for calculation of block borders...might change). Defaults to chrom-sizes for hg19 from UCSC.") #TODO: maybe use VinDenAlt manifesto here; later:
	parser.add_argument("--centro", help="optional: bed-file containing chr,start,stop for the centromere of the chromosome. If defined the centromere will be excluded.")

	args = parser.parse_args()

	## vcf from stdin
	#if sys.stdin.isatty():
		#sys.exit("Error: Needs vcf from stdin, usage: 'zcat file.vcf.gz | d_stats.py pop1 pop2 pop3 pop4 -i pop-info -n 20 -t'")
	#vcf = sys.stdin

	########

	## get population matches pop1, pop2, pop3, pop4
	pw_pops = D.get_pops_from_files(args.pop1, args.pop2, args.pop3, args.pop4)
	#D.print_table(pw_pops)
	#D.print_table_to_file(pw_pops, "pw_pops") # TODO: is that needed? was active in original d_stats.py

	## TODO: may get centro range and chrom range when first line is parsed (chrom-arg not needed anymore)
	## get centromere range
	if args.centro:
		centro_range = D.get_range(args.chrom ,args.centro)
	else:
	    centro_range = [0,0]
	print(centro_range)

	## NEW: avoid chrom-ranges and block borders - only take x-kb steps when going through vcf
	### get chrom range for blocking if nr of blocks is defined, and not blocksize
	#chrom_range = D.get_range(args.chrom ,args.startstop)
	#print(chrom_range)
	## get borders for chromosome-blocking
	#block_borders = D.get_block_borders(chrom_start, chrom_end, centro_range[0], centro_range[1], n_blocks, kb)
	#with open("block_borders","w") as blockout:
		#blockout.write("\n".join(map(str,block_borders)))




## get population and sex for every vcf-individual
#pop_dict = F.get_pop_dict(sample_info)


## get vcf-columns for every population (NEW: and for African2)
#vcf_header = vcf_input.readline().rstrip().split()
##print vcf_header
#pop_colnums, col_gender = F.get_pop_colnumbers(pop_dict, vcf_header)
##print pop_colnums
#print "Number of individuals in Africa2: ", len(pop_colnums["Africa2"])



## compute blockwise ABBA  [ [[abba][baba]] [[abba][baba]] ...] (NEW param-Afr fixed)
#blocks, sites_comp = F.abba_block_sums(vcf_input, pop_colnums, pw_pops, col_gender, block_borders, centro_range, private, transver, affixed)

### rearrange output and print to file [[pop1][pop2][pop3][block][abba][baba]]
#blocks_out = F.rearrange_blocks(pw_pops, blocks)
#with open("out_blocks","w") as out:
	#out.write('\t'.join(["pop1","pop2","pop3","pop4","block","abba","baba"])+"\n")
#F.print_table_to_file(blocks_out, "out_blocks", mode="a")

### NEW print sites per comparison to file [[pop1][pop2][pop3][sites]]
#pw_pops.append(sites_comp)
#with open("sites_comp","w") as out:
	#out.write('\t'.join(["pop1","pop2","pop3","pop4","n_sites"])+"\n")
#F.print_table_to_file(pw_pops, "sites_comp", mode="a")




if __name__ == "__main__":
	main()



