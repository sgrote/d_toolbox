#!/usr/bin/python3

''' 
Filter one vcf file with header, from STDIN
filtering for lines and individual genotypes with condition-string.
Header and QUAL, FILTER, INFO are kept, FORMAT is kept or "GT" if gt_only.
Writes to STDOUT


generate test-file
cd ~/Test
CTEAM=/mnt/sequencedb/gendivdata/2_genotypes/human/SGDP/SGDP_v3_May2016/combined_vcf/c_team_chr21.vcf.gz
zcat $CTEAM | head -n 10000 | cut -f -50 | awk '$6<1 {$7="LowQual"}; $6>100 {$7="PASS"}  1' OFS="\t" > test_lines_SGDP

'''

import argparse
import sys

import vcf_line_filter as F


def main():
	parser = argparse.ArgumentParser(description='Filter one vcf file, from STDIN filtering for lines and individual genotypes with condition-string. Header and QUAL, FILTER, INFO are kept, FORMAT is kept or "GT" if gt_only. Writes to STDOUT.', usage='zcat file1.vcf.gz | vcf_filter.py --filter "QUAL > 3 and AC >= 10" --filter_ind "GF >= 1"')
	parser.add_argument("--filter", help="Expression to filter vcf1, e.g. 'QUAL > 0 and AC > 10'. Acts on QUAL, FILTER and INFO fields. CAUTION: needs whitespaces around keywords.")
	parser.add_argument("--filter_ind", help="Expression to filter individual genotypes of vcf1, e.g. 'GF >= 0 or GT == 1/1'. Acts on keywords in FORMAT. Not passing genotypes will be replaced by './.' . CAUTION: needs whitespaces around keywords.")
	parser.add_argument("--var", action="store_true", help="Restrict to variable sites")
	parser.add_argument("--gt_only", action="store_true", help="Restrict to genotypes in sample columns")
	parser.add_argument("--keep_miss", action="store_true", help="Keep lines that were filtered out as missing data ./.")
	
	args = parser.parse_args()

	if sys.stdin.isatty():
		sys.exit('Error: Needs vcf from stdin, e.g. "zcat file1.vcf.gz | vcf_filter.py --filter "QUAL > 3 and AC >= 10" --filter_ind "GF >= 1""')
	vcf = sys.stdin

	
	### 1) Header
	
	# vcf header, keep line with donor-ids
	line = vcf.readline()
	while line[:6] != "#CHROM":
		sys.stdout.write(line)
		line = vcf.readline()
	# write sample-line
	sys.stdout.write(line) 


	#### 2) Genotypes
		
	for line in vcf:
	
		v = line.rstrip().split()
		
		## for every line: check if line passes and filter on individual genotypes
		try:
			# returns (modified) line if it passes filter, if not None or all ./. line (if keep_miss)
			out = F.combi_filter(v, filter1=args.filter, filter_ind1=args.filter_ind, gt_only=args.gt_only, var=args.var, keep_miss=args.keep_miss)
			if out:
				sys.stdout.write("\t".join(out) + "\n")
		except (IndexError, ValueError) as errore:
			print(errore)
			print(v)


if __name__ == "__main__":
	main()


''' test
FILTER=/mnt/expressions/steffi/D/d_toolbox/vcf_filter.py
cat test_lines_SGDP | $FILTER --filter "QUAL > 0" | less -S
cat test_lines_SGDP | $FILTER --filter "QUAL > 0" --keep_miss | less -S
cat test_lines_SGDP | $FILTER --filter "QUAL > 0" --filter_ind "GF >= 0" | less -S
cat test_lines_SGDP | $FILTER --filter_ind "GF >= 0" | less -S
# keep only lines where ALT allele is present after filtering individual genotypes
cat test_lines_SGDP | $FILTER --filter_ind "GF >= 0" --var | less -S

'''




