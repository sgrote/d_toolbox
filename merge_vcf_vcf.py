#!/usr/bin/python3

'''
write to STDOUT for multiple merging steps
take header from one of the files or add a minimal standard header (default)
QUAL, FILTER, INFO fields are not merged/updated but replaced by "."
FORMAT is taken from one of the files, or if gt_only replaced by "GT"
optional filtering on quality (for an already merged file it's harder to filter on file-specific criteria)
'''

# maybe TODO: accurate merging of FORMAT, "1/1:.:." for missing (but mostly --gt_only --> not needed)
# maybe TODO: modify merge_alt_alleles to take position of genotype as input (now for "1/1:other:info the genotype is assumed to be the first entry")


import sys
import argparse
import gzip

import merge_vcf_GTs as M
import vcf_line_filter as F


### helper

# If only one line is present or passes line-filter
def filter_and_fill(line1, n_donors2, filter1=None, filter_ind1=None, gt_only=False, fill2=False, var=None, keep_miss=None):
	''' 
	Filter line1 with condition string, if pass also filter individual genotypes and create ./. or 0/0 genotypes for line2. Also remove QUAL, FILTER and INFO
	Out: (filtered_line1, filled_gt2) or (None, None) if filter not passed
	'''
	# skip invariable sites even when keep_miss
	if var and line1[4] == ".":
		return None, None
	# create genotypes for line2:
	if fill2:
		gt2 = ['0/0'] * n_donors2
	else:
		gt2 = ['./.'] * n_donors2
	# filter line1
	if filter1 and not F.filter_line(filter1, line1):
		if keep_miss:
			out_line = line1[:5] + [".",".",".","GT"] + ["./."] * (len(line1)-9)
			return out_line, gt2
		else:
			return None, None
	# filter line1 genotypes
	out_line = F.filter_gt(line1, filter_ind1, gt_only)
	# check if any of the (ALT)-alleles is left after gt-filterig
	if not keep_miss and not F.check_gt(out_line, var):
		return None, None
	# remove QUAL FILTER INFO for merged line
	out_line[5:8] = [".",".","."]
	return out_line, gt2

''' test
line1 = ['21','148','.','C','T','10','LowQual', 'AC=27;AN=558;AF=a_string','GT:GF','0/1:-1','0/0:3','./0:-1']
line_clean = ['21','148','.','C','T','.','.','.','GT:GF','0/1:-1','0/0:3','./0:-1']
filter_and_fill(line1, 2) == (line_clean, ['./.', './.'])
filter_and_fill(line1, 2, "QUAL > 9") == (line_clean, ['./.', './.'])
filter_and_fill(line1, 3, "QUAL > 9", fill2=True) == (line_clean, ['0/0', '0/0', '0/0'])
filter_and_fill(line1, 3, "QUAL > 10", fill2=True) == (None, None)
filter_and_fill(line1, 2, "QUAL > 9 and FILTER == 'LowQual'") == (line_clean, ['./.', './.'])
filter_and_fill(line1, 2, "QUAL > 9 and FILTER != 'LowQual'") == (None, None)
filter_and_fill(line1, 2, "AN > 500 or AC > 30") == (line_clean, ['./.', './.'])
filter_and_fill(line1, 2, "AN > 500 and AC > 30") == (None, None)
filter_and_fill(line1, 1, "AN > 500", gt_only=True) == (line_clean[:8]+['GT','0/1','0/0','./0'], ['./.'])
filter_and_fill(line1, 2, "QUAL > 0", "GF >= 0") == (line_clean[:9]+['./.:-1','0/0:3','./.:-1'], ['./.', './.'])
filter_and_fill(line1, 1, "QUAL > 0", "GF >= 0", gt_only=True) == (line_clean[:8] + ['GT','./.','0/0','./.'],['./.'])
filter_and_fill(line1, 1, "QUAL > 0", "GF >= 0", gt_only=True, var=True) == (None, None)
filter_and_fill(line1, 1, "QUAL > 0", "GF >= 0", gt_only=True, var=True, keep_miss=True) == (line_clean[:8] + ['GT','./.','0/0','./.'],['./.'])
filter_and_fill(line1, 1, filter_ind1="GF > 3", gt_only=True) == (None, None)
filter_and_fill(line1, 1, filter_ind1="GF > 3", gt_only=True, keep_miss=True) == (line_clean[:8]+['GT','./.','./.','./.'], ['./.'])
filter_and_fill(line1, 1, filter_ind1="GF > 3", gt_only=True, keep_miss=True, fill2=True) == (line_clean[:8]+['GT','./.','./.','./.'], ['0/0'])
filter_and_fill(['21','148','.','C','.','10','.', '.','GT:GF','0/0:-1'],2,var=True) == (None, None)
filter_and_fill(['21','148','.','C','.','10','.', '.','GT:GF','0/0:-1'],2,var=True, keep_miss=True) == (None, None)

'''


# If both lines are present and pass line-filter
def merge_and_filter(line1, line2, filter_ind1=None, filter_ind2=None, gt_only=False, var=False, keep_miss=False):
	''' 
	When positions match and both lines pass line-filter:
	merge, filter genotypes, check var
	'''
	# avoid merging if var not passed
	if var and line1[4] == "." and line2[4] == ".":
		return None
	line1 = F.filter_gt(line1, filter_ind1, gt_only)
	line2 = F.filter_gt(line2, filter_ind2, gt_only)
	out = M.merge_lines(line1, line2)
	# check if line is still valid after filtering individual genotypes
	if not keep_miss and not F.check_gt(out, var):
		return None
	else:
		return out

''' test
line1 = ['21','148','.','C','T','10','.','AC=27','GT:GF','0/1:3','./0:-1']
line2 = ['21','148','.','C','A','10','.','.','GT','0|1','0|1']
line3 = ['21','148','.','C','.','10','PASS','.','GT','0|0','0|0']
merge_and_filter(line1, line3) == ['21','148','.','C','T','.','.','.','GT:GF','0/1:3','./0:-1','0|0','0|0']
merge_and_filter(line3, line3, var=True) == None
merge_and_filter(line3, line3, var=True, keep_miss=True) == None
merge_and_filter(line1, line2) == ['21','148','.','C','A,T','.','.','.','GT:GF','0/2:3','./0:-1','0|1','0|1']
merge_and_filter(line1, line3, filter_ind1="GF > 3") == ['21','148','.','C','T','.','.','.','GT:GF','./.:3','./.:-1','0|0','0|0']
merge_and_filter(line1, line3, filter_ind1="GF > 3", var=True) == None
merge_and_filter(line1, line3, filter_ind1="GF > 3", var=True, keep_miss=True, gt_only=True) == ['21','148','.','C','T','.','.','.','GT','./.','./.','0|0','0|0']

'''

def main():
	parser = argparse.ArgumentParser(description='Merge two vcf files, one from STDIN the other as first command-line argument (.vcf or .vcf.gz). Specify whether missing positions in one file should be converted to "./." or "0/0" with fill_ref1, fill_ref2 flags. Optional filtering for lines and individual genotypes with condition-string. Header is either a minimal default header or taken from one of the files. QUAL, FILTER, INFO are replaced by ".", FORMAT is taken from one of the files or "GT" if gt_only. Writes to STDOUT.', usage='zcat file1.vcf.gz | merge_vcf_vcf.py file2.vcf.gz --fill1 --filter2 "QUAL > 3 and AC >= 10"')
	parser.add_argument("vcf_2_file", help="Second vcf file to merge (.vcf or .vcf.gz)")
	parser.add_argument("--var", action="store_true", help="Restrict to variable sites")
	parser.add_argument("--gt_only", action="store_true", help="Restrict to genotypes in sample columns")
	parser.add_argument("--fill1", action="store_true", help="Fill vcf1 with REF 0/0 for missing positions, else ./.")
	parser.add_argument("--fill2", action="store_true", help="Like --fill1 for vcf2")
	parser.add_argument("--require1", action="store_true", help="Restrict to positions present in vcf1")
	parser.add_argument("--require2", action="store_true", help="Restrict to positions present in vcf2")
	parser.add_argument("--keep_miss", action="store_true", help="Keep lines that were filtered out as missing data ./.")
	parser.add_argument("--filter1", help="Expression to filter vcf1, e.g. 'QUAL > 0 and AC > 10'. Acts on QUAL, FILTER and INFO fields. CAUTION: needs whitespaces around keywords.")
	parser.add_argument("--filter2", help="Like --filter1 for vcf2")
	parser.add_argument("--filter_ind1", help="Expression to filter individual genotypes of vcf1, e.g. 'GF >= 0 or GT == 1/1'. Acts on keywords in FORMAT. Not passing genotypes will be replaced by './.' . CAUTION: needs whitespaces around keywords.")
	parser.add_argument("--filter_ind2", help="Like --filter_ind1 for vcf2")
	parser.add_argument("--header", choices=["vcf1", "vcf2"], help="optional: take the header from vcf1 or vcf2. If not specified a minimal standard header is added.")
	
	args = parser.parse_args()

	if sys.stdin.isatty():
		sys.exit("Error: Needs vcf1 from stdin, e.g. 'zcat file1.vcf.gz | merge_vcf_vcf.py file2.vcf'")
	vcf_1 = sys.stdin

	if args.vcf_2_file.endswith('.gz'):
		opener = gzip.open
	else:
		opener = open
	with opener(args.vcf_2_file, 'rt') as vcf_2:	
		
		### 1) Header
		
		# vcf_2 header, keep line with donor-ids
		v2 = vcf_2.readline()
		while v2[:6] != "#CHROM":
			if args.header=="vcf2":
				sys.stdout.write(v2)
			v2 = vcf_2.readline()
		v2_donors = v2.rstrip().split()[9:]	
		
		# vcf_1 header, keep line with donor-ids
		v1 = vcf_1.readline()
		while v1[:6] != "#CHROM":
			if args.header=="vcf1":
				sys.stdout.write(v1)
			v1 = vcf_1.readline()

		# add default header
		if not args.header:
			sys.stdout.write('##fileformat=VCFv4.1' + "\n")
			sys.stdout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype calls or randomly sampled reads (for low-coverage archaics)">' + "\n")
			sys.stdout.write('##reference=file:/mnt/solexa/Genomes/hg19_evan/whole_genome.fa' + "\n")
		
		# print combined donors
		sys.stdout.write("\t".join([v1.rstrip()] + v2_donors) + "\n")

		#### 2) Genotypes
			
		v1 = vcf_1.readline().split()
		v2 = vcf_2.readline().split()
		if v1[0] != v2[0]:
			sys.exit("Error: Chromosomes do not match: vcf_1-chrom={}, vcf_2-chrom={}.".format(v1[0], v2[0]))
		## number of donors
		n_v1 = len(v1)-9
		n_v2 = len(v2)-9
		### iterate over lines from both files and merge, until end of both files is reached
		while not((len(v1) == 0) and (len(v2) == 0)):
			## for every line: check if line passes and filter on individual genotypes
			try:
				out = None
				# pos not in v2, including the case that v2 file has reached end
				# -> if v2 not required and v1 passes filter: add v2 0/0 or ./. GTs and read next v1 line
				if (len(v2) == 0) or (len(v1) != 0 and int(v1[1]) < int(v2[1])):
					if not args.require2:
						v1, gt2 = filter_and_fill(v1, n_v2, args.filter1, args.filter_ind1, args.gt_only, args.fill2, args.var, args.keep_miss)
						if v1:
							out = v1 + gt2
					v1 = vcf_1.readline().split()
				## pos not v1
				elif (len(v1) == 0) or (len(v2) != 0 and int(v2[1]) < int(v1[1])):
					if not args.require1:
						v2, gt1 = filter_and_fill(v2, n_v1, args.filter2, args.filter_ind2, args.gt_only, args.fill1, args.var, args.keep_miss)
						if v2:
							out = v2[:9] + gt1 + v2[9:]
					v2 = vcf_2.readline().split()
				## positions are present in both vcfs	
				elif v1[1] == v2[1]:
					# (NEW: also check that REF-alleles match (no-match seen in indels in 1000 genomes))
					if v1[3] == v2[3]:
						v1_pass = F.check_filter(v1, args.filter1)
						v2_pass = F.check_filter(v2, args.filter2)
						if v1_pass and not v2_pass:
							# fill2=False because v2 did not pass filter and should be ./. and not 0/0
							v1, gt2 = filter_and_fill(v1, n_v2, None, args.filter_ind1, args.gt_only, False, args.var, args.keep_miss)
							if v1: # could fail despite v1_pass when individual GTs don't pass or invariable
								out = v1 + gt2
						elif v2_pass and not v1_pass:
							v2, gt1 = filter_and_fill(v2, n_v1, None, args.filter_ind2, args.gt_only, False, args.var, args.keep_miss)
							if v2:
								out = v2[:9] + gt1 + v2[9:]
						# both pass line-filter: merge, filter individual genotypes and check var again
						elif v1_pass and v2_pass:
							out = merge_and_filter(v1, v2, args.filter_ind1, args.filter_ind2, args.gt_only, args.var)
					# read next lines from both files (also if not (v1_pass or v2_pass))
					v1 = vcf_1.readline().split()
					v2 = vcf_2.readline().split()
				# print merged line if any
				if out:
					sys.stdout.write("\t".join(out) + "\n")
			except (IndexError, ValueError) as errore:
				print(errore)
				print(v1)
				print(v2)


if __name__ == "__main__":
	main()


''' test
VCF1=/mnt/scratch/steffi/D/Vcfs/mergedArchModernApes/merged_high_chr21.vcf.gz
VCF2=/mnt/sequencedb/gendivdata/2_genotypes/human/SGDP/SGDP_v3_May2016/combined_vcf/c_team_chr21.vcf.gz
MERGE=/mnt/expressions/steffi/D/d_toolbox/merge_vcf_vcf.py
zcat $VCF1 | $MERGE $VCF2 | less -S
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter2 "QUAL > 0" | less -S
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter2 "QUAL > 0" --filter_ind2 "GF >= 1" | less -S
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter2 "QUAL > 0" --filter_ind2 "GF >= 1" --keep_miss | less -S
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter_ind2 "GF >= 1" | less -S
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter_ind2 "GF >= 1" --keep_miss | less -S
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter2 "QUAL > 0" --gt_only | less -S
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter1 "QUAL > 70" --gt_only | less -S
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter2 "QUAL > 0" --gt_only --var | less -S
zcat $VCF1 | $MERGE $VCF2 --filter_ind1 "GT == '1/1'" --filter2 "QUAL > 0" | less -S
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter2 "QUAL > 0" --filter_ind2 "GF >= 1" --gt_only --var | less -S
zcat $VCF1 | $MERGE $VCF2 --gt_only --var --require1 | less -S
zcat $VCF1 | $MERGE $VCF2 --gt_only --var --require2 | less -S
zcat $VCF1 | $MERGE $VCF2 --gt_only --var --require1 --require2 | less -S

# require archaics to be present (since they are already filtered and needed in all following D-stats)
echo `date`
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter2 "QUAL > 0" --filter_ind2 "GF >= 1" --gt_only --var --require1 > test_merge4.vcf
echo `date`

'''



