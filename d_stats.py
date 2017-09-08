#!/usr/bin/python

'''
Calculate blockwise ABBA and BABA and write to file "out_blocks"
containing lines with the 3 pops and ABBA, BABA for every block
NEW: read pop1-pop4 from file, don't use chimp as outgroup anymore
'''


import sys
import imp

print sys.argv

# variable input
sample_info = sys.argv[1]	# donor metadata
chrom = sys.argv[2]			# chromosome
chrom_start = sys.argv[3]	# start informative sites
chrom_end = sys.argv[4]		# end informative sites
n_blocks = sys.argv[5]		# Nr. of blocks
kb = sys.argv[6]			# should n_blocks be taken as length in kb rather then nr of blocks?
private = sys.argv[7]		# should Denisovan be required to be homozygous ancestral?
transver = sys.argv[8]		# should only transversions be included?
affixed = sys.argv[9]		# fixed in Africa2?
pop1 = sys.argv[10]			# file with pop1s in rows (as in vcf header)
pop2 = sys.argv[11]
pop3 = sys.argv[12]
pop4 = sys.argv[13]

vcf_input = sys.stdin		# merged archaic and modern 'genotypes'


F = imp.load_source('d_functions', '/mnt/expressions/steffi/D/scripts_d_mateja/d_functions.py')


# get centromere range
#centro_range = F.get_range(chrom ,"/mnt/expressions/steffi/centromereRanges/centromere_ranges")
centro_range = [0,0]
# get borders for chromosome-blocking
block_borders = F.get_block_borders(chrom_start, chrom_end, centro_range[0], centro_range[1], n_blocks, kb)
with open("block_borders","w") as blockout:
	blockout.write("\n".join(map(str,block_borders)))
# get population and sex for every vcf-individual
pop_dict = F.get_pop_dict(sample_info)

# get vcf-columns for every population (NEW: and for African2)
vcf_header = vcf_input.readline().rstrip().split()
#print vcf_header
pop_colnums, col_gender = F.get_pop_colnumbers(pop_dict, vcf_header)
#print pop_colnums
print "Number of individuals in Africa2: ", len(pop_colnums["Africa2"])

# get population matches pop1, pop2, pop3, pop4
pw_pops = F.get_pops_from_files(pop1, pop2, pop3, pop4)
#F.print_table(pw_pops)
F.print_table_to_file(pw_pops, "pw_pops")

# compute blockwise ABBA  [ [[abba][baba]] [[abba][baba]] ...] (NEW param-Afr fixed)
blocks, sites_comp = F.abba_block_sums(vcf_input, pop_colnums, pw_pops, col_gender, block_borders, centro_range, private, transver, affixed)

## rearrange output and print to file [[pop1][pop2][pop3][block][abba][baba]]
blocks_out = F.rearrange_blocks(pw_pops, blocks)
with open("out_blocks","w") as out:
	out.write('\t'.join(["pop1","pop2","pop3","pop4","block","abba","baba"])+"\n")
F.print_table_to_file(blocks_out, "out_blocks", mode="a")

## NEW print sites per comparison to file [[pop1][pop2][pop3][sites]]
pw_pops.append(sites_comp)
with open("sites_comp","w") as out:
	out.write('\t'.join(["pop1","pop2","pop3","pop4","n_sites"])+"\n")
F.print_table_to_file(pw_pops, "sites_comp", mode="a")





