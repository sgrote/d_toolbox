#!/usr/bin/python3

'''
write to STDOUT for multiple merging steps
header is not merged but only last line with sample names is kept
QUAL, FILTER, INFO fields are not merged/updated but replaced by "."
optional filtering on quality (for an already merged file it's harder to filter on file-specific criteria)
'''

# TODO: accurate merging of FORMAT, "1/1:.:." for missing (but mostly --gt_only --> not needed)
# maybe TODO: merge headers
# maybe TODO: modify merge_alt_alleles to take position of genotype as input (now for "1/1:other:info the genotype is assumed to be the first entry")


import sys
import imp
import argparse
import gzip

import merge_vcf_GTs as M
import vcf_line_filter as F


# helper
# TODO: check if any o the ALT-alleles is left after gt-filterig; if --var remove site if not
def filter_and_fill(line1, n_donors2, filter1=None, filter_ind1=None, gt_only=False, fill2=False):
	''' 
	Filter line1 with condition string, if pass also filter individual genotypes and create ./. or 0/0 genotypes for line2. Also remove QUAL, FILTER and INFO
	Out: (filtered_line1, filled_gt2) or (None, None) if filter not passed
	'''
	# filter line1
	if filter1 and not F.filter_line(filter1, line1):
		return None, None
	# filter line1 genotypes
	out_line = F.filter_gt(line1, filter_ind1, gt_only)
	# remove QUAL FILTER INFO for merged file
	out_line[5:8] = [".",".","."]
	# create genotypes for line2:
	if fill2:
		gt2 = ['0/0'] * n_donors2
	else:
		gt2 = ['./.'] * n_donors2
	return out_line, gt2

''' test
line1 = ['21','148','.','C','T','10','LowQual', 'AC=27;AN=558;AF=a_string','GT:GF','0/1:-1','0/1:3','./0:-1']
line_clean = ['21','148','.','C','T','.','.','.','GT:GF','0/1:-1','0/1:3','./0:-1']
filter_and_fill(line1, 2) == (line_clean, ['./.', './.'])
filter_and_fill(line1, 2, "QUAL > 9") == (line_clean, ['./.', './.'])
filter_and_fill(line1, 3, "QUAL > 9", fill2=True) == (line_clean, ['0/0', '0/0', '0/0'])
filter_and_fill(line1, 3, "QUAL > 10", fill2=True) == (None, None)
filter_and_fill(line1, 2, "QUAL > 9 and FILTER == 'LowQual'") == (line_clean, ['./.', './.'])
filter_and_fill(line1, 2, "QUAL > 9 and FILTER != 'LowQual'") == (None, None)
filter_and_fill(line1, 2, "AN > 500 or AC > 30") == (line_clean, ['./.', './.'])
filter_and_fill(line1, 2, "AN > 500 and AC > 30") == (None, None)
filter_and_fill(line1, 1, "AN > 500", gt_only=True) == (line_clean[:8]+['GT','0/1','0/1','./0'], ['./.'])
filter_and_fill(line1, 2, "QUAL > 0", "GF >= 0") == (line_clean[:9]+['./.:-1','0/1:3','./.:-1'], ['./.', './.'])
filter_and_fill(line1, 1, "QUAL > 0", "GF >= 0", gt_only=True) == (line_clean[:8] + ['GT','./.','0/1','./.'],['./.'])

'''



def main():
	parser = argparse.ArgumentParser(description='Merge two vcf files, one from STDIN the other as first command-line argument (.vcf or .vcf.gz). Specify whether missing positions in one file should be converted to "./." or "0/0" with fill_ref1, fill_ref2 flags. Optional filtering for lines and individual genotypes with condition-string. Header is reduced to sample names. Writes to STDOUT.', usage='zcat file1.vcf.gz | merge_vcf_vcf.py file2.vcf.gz --fill1 --filter2 "QUAL > 3 and AC >= 10"')
	parser.add_argument("vcf_2_file", help="Second vcf file to merge (.vcf or .vcf.gz)")
	parser.add_argument("--var", action="store_true", help="Restrict to variable sites")
	parser.add_argument("--gt_only", action="store_true", help="Restrict to genotypes in sample columns")
	parser.add_argument("--fill1", action="store_true", help="Fill vcf1 with REF 0/0 for missing positions, else ./.")
	parser.add_argument("--fill2", action="store_true", help="Like --fill1 for vcf2")
	parser.add_argument("--filter1", help="Expression to filter vcf1, e.g. 'QUAL > 0 and AC > 10'. Acts on QUAL, FILTER and INFO fields. CAUTION: needs whitespaces around keywords.")
	parser.add_argument("--filter2", help="Like --filter1 for vcf2")
	parser.add_argument("--filter_ind1", help="Expression to filter individual genotypes of vcf1, e.g. 'GF >= 0 or GT == 1/1'. Acts on keywords in FORMAT. Not passing genotypes will be replaced by './.' . CAUTION: needs whitespaces around keywords.")
	parser.add_argument("--filter_ind2", help="Like --filter_ind1 for vcf2")
	
	args = parser.parse_args()

	if sys.stdin.isatty():
		print("Needs vcf1 from stdin, e.g. 'zcat file1.vcf.gz | merge_vcf_vcf.py file2.vcf'")
	vcf_1 = sys.stdin

	if args.vcf_2_file.endswith('.gz'):
	    opener = gzip.open
	else:
	    opener = open
	with opener(args.vcf_2_file, 'rt') as vcf_2:	
		
		### 1) Header
		
		# skip vcf_2 header, only keep line with donor-ids
		v2 = vcf_2.readline()
		while v2[:6] != "#CHROM":
			v2 = vcf_2.readline()
		v2_donors = v2.rstrip().split()[9:]	
		
		# skip vcf_1 header, only keep line with donor-ids
		v1 = vcf_1.readline()
		while v1[:6] != "#CHROM":
			v1 = vcf_1.readline()
		sys.stdout.write("\t".join([v1.rstrip()] + v2_donors) + "\n")
	
		#### 2) Genotypes
			
		v1 = vcf_1.readline().split()
		v2 = vcf_2.readline().split()
		## number of donors
		n_v1 = len(v1)-9
		n_v2 = len(v2)-9
		### iterate over lines from both files and merge, until end of both files is reached
		while not((len(v1) == 0) and (len(v2) == 0)):
			## for every line: check if line passes and filter on individual genotypes
			try:
				out = None
				# pos not in v2, including the case that v2 file has reached end
				# -> if v1 passes filter add v2 0|0 GTs to v1 and read next v1 line
				if (len(v2) == 0) or (len(v1) != 0 and int(v1[1]) < int(v2[1])):
					v1, gt2 = filter_and_fill(v1, n_v2, args.filter1, args.filter_ind1, args.gt_only, args.fill2)
					if v1:
						out = v1 + gt2
					v1 = vcf_1.readline().split()
				## pos not v1
				elif (len(v1) == 0) or (len(v2) != 0 and int(v2[1]) < int(v1[1])):
					v2, gt1 = filter_and_fill(v2, n_v1, args.filter2, args.filter_ind2, args.gt_only, args.fill1)
					if v2:
						out = v2[:9] + gt1 + v2[9:]
					v2 = vcf_2.readline().split()
				## positions are present in both vcfs	
				elif v1[1] == v2[1]:
					v1_pass = F.check_filter(v1, args.filter1)
					v2_pass = F.check_filter(v2, args.filter2)
					if v1_pass and not v2_pass:
						v1, gt2 = filter_and_fill(v1, n_v2, None, args.filter_ind1, args.gt_only, args.fill2)
						out = v1 + gt2
					elif v2_pass and not v1_pass:
						v2, gt1 = filter_and_fill(v2, n_v1, None, args.filter_ind2, args.gt_only, args.fill1)
						out = v2[:9] + gt1 + v2[9:]
					elif v1_pass and v2_pass:
						v1 = F.filter_gt(v1, args.filter_ind1, args.gt_only)
						v2 = F.filter_gt(v2, args.filter_ind2, args.gt_only)
						out = M.merge_lines(v1, v2)	
					# read next lines from both files (also if not (v1_pass or v2_pass))
					v1 = vcf_1.readline().split()
					v2 = vcf_2.readline().split()
				# print merged line if any
				if out and not (args.var and out[4]=="."):
					sys.stdout.write("\t".join(out) + "\n")
			except (IndexError, ValueError) as errore:
				print(errore)
				print(v1)
				print(v2)


if __name__ == "__main__":
	main()


''' test
VCF1=/mnt/scratch/steffi/D/Vcfs/Altai_Vindija_Denis/altai_vindija_denis_chr21_filtered_apes.vcf.gz
VCF2=/mnt/sequencedb/gendivdata/2_genotypes/human/SGDP/SGDP_v3_May2016/combined_vcf/c_team_chr21.vcf.gz
MERGE=/mnt/expressions/steffi/D/d_toolbox/merge_vcf_vcf.py
zcat $VCF1 | $MERGE $VCF2 | less -S
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter2 "QUAL > 0" | less -S
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter2 "QUAL > 0" --filter_ind2 "GF >= 1" | less -S
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter2 "QUAL > 0" --gt_only | less -S
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter1 "QUAL > 70" --gt_only | less -S
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter2 "QUAL > 0" --gt_only --var | less -S
zcat $VCF1 | $MERGE $VCF2 --filter_ind1 "GT == '1/1'" --filter2 "QUAL > 0" | less -S
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter2 "QUAL > 0" --filter_ind2 "GF >= 1" --gt_only --var | less -S

echo `date`
zcat $VCF1 | $MERGE $VCF2 --fill2 --filter_ind2 "GF >= 1" --gt_only --var > test_merge.vcf
echo `date`

'''


