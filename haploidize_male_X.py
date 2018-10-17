#!/usr/bin/python3

'''
go through a X-chrom VCF and
make all male genotypes outside pseudo autosomal region haploid

get sample info from file
take vcf from stdin
PAR coords are hard-coded

'''

import sys
import argparse
import pandas as pd

import d_functions as D


def main():
	parser = argparse.ArgumentParser(description='Takes pop-info from file, vcf from stdin and replaces diploid male GTs outside PAR-region by haploid GT. Hets are masked as ".". Writes to Stdout', usage='zcat chrX.vcf.gz | haploidize_male_X.py -i pop_info.csv')
	# mandatory
	parser.add_argument("-i", "--info", required=True, help="mandatory: csv-file with population metadata for pops. Needs columns [sample; population; sex]; with sample being identical to the name in the vcf.")

	args = parser.parse_args()

	# vcf from stdin
	if sys.stdin.isatty():
		sys.exit("Error: Needs vcf from stdin, usage: 'zcat chrX.vcf.gz | haploidize_male_X.py -i pop_info.csv")
	vcf = sys.stdin

	########

	# vcf header, keep line with donor-ids
	line = vcf.readline()
	while line[:6] != "#CHROM":
		sys.stdout.write(line)
		line = vcf.readline()
	# write sample-line
	sys.stdout.write(line) 

	# from info-file and vcf-header: getsex for every col
	male_cols = get_sex_cols(args.info, line.rstrip().split())
	
	# Replace genotypes in following lines
	haploidize(vcf, male_cols)



def get_sex_cols(info_file, vcf_header, sex="male"):
	'''
	in: info_file csv-file with columns [sample, pop, sex]
		vcf_header as list
	out: list [colnumbers of male samples]
	'''
	male_cols = []
	info = pd.read_csv(info_file)
	# for each donor-id from vcf-header, get sex if present in info-file
	for i in range(9, len(vcf_header)):
		if any(info['sample'] == vcf_header[i]):
			sample_info = info[info['sample'] == vcf_header[i]]
			# get sex of col-number 
			# (first entry of sample_info['sex'], since this column is all the same for one sample)
			# (one sample can have multiple entires in info since it belongs to different populations)
			sample_sex = sample_info['sex'].iloc[0]
			if sample_sex == sex:
				male_cols.append(i)
	return male_cols

''' test
vcf_header = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','S_Esan-1','S_Esan-2','S_Luo-1']
info_file = '/mnt/expressions/steffi/D/infofiles/combined_info.csv'
male_cols = get_sex_cols(info_file, vcf_header)
male_cols == [9, 11]
female_cols = get_sex_cols(info_file, vcf_header, 'female')
female_cols == [10]

'''


def haploidize_line(vcf_line, male_cols):
	'''
	haploidize genotypes of male_cols in vcf_line
	'''
	for colnumber in male_cols:
		gt = vcf_line[colnumber]
		# leave haploid GTs
		if len(gt) == 1:
			continue
		# homozygous diploids -> haploid;  heterozygous diploids -> mask  
		elif len(gt) == 3:
			alleles = set([gt[i] for i in [0,2]])
			if len(alleles) == 1:
				vcf_line[colnumber] = alleles.pop()
			else:
				vcf_line[colnumber] = "."
		else:
			sys.stderr.write("unexpected genotype:\n")
			sys.stderr.write(gt)
			#sys.exit()
	return vcf_line


''' test
line = ["0/0","0/1","1","2","2|2","1/1","1/2","0|0"]
cols = [0,2,4,5,6,7]
haploidize_line(line, cols) == ["0","0/1","1","2","2","1",".","0"]

'''



def haploidize(vcf, male_cols):
	'''
	haploidize genotypes of male_cols in vcf
	CAUTION this assumes pure genotype-input
	'''
	for line in vcf:
		try:
			line = line.rstrip().split()
			pos = int(line[1])
			# leave PAR as it is
			if (pos > 60000 and pos < 2699520) or (pos > 154931043 and pos < 155260560):
				sys.stdout.write("\t".join(line) + "\n")
			# loop over all male samples in vcf-file and add genotype to dictionary
			else:
				haplo_line = haploidize_line(line, male_cols)
				sys.stdout.write("\t".join(haplo_line) + "\n")
			
		except (IndexError, ValueError) as errore:
			print(errore)
			print(line)
	

if __name__ == "__main__":
	main()





