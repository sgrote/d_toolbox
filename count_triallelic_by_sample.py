#!/usr/bin/python3

'''
Given a file with samples
count sites in VCF that are triallelic due to those samples only (and biallelic without them)
write a table with
(a) counts per individual (singletons)
(b) all triallelics caused by one or more of the samples (including all (a))
(c) triallic counts outside input samples
(d) biallelic counts outside input samples (including all (b), might be biallelic in samples, too)

(biallelics caused only by input samples are not counted anywhere)

Purpose: get an idea of how many true biallelic sites are triallelic by
(potential) errors in low-coverage genomes or long-branch ape reference genomes

write to file
take vcf from stdin
'''

import sys
import argparse

import d_functions as D


def main():
	parser = argparse.ArgumentParser(description='Given a file with samples, count sites in VCF that are triallelic due to those samples only. Write a table with counts per individual (singletons) including a line with counts where triallelic is caused by more than one of the samples (but still all non-input samples would be biallelic). For comparison also add to lines with counts of biallelic and triallelic sites outside the samples.', usage='zcat chr21.vcf.gz | count_triallelic_by_sample.py sample_file')
	parser.add_argument("samples", type=str, help="File with sample names like in VCF header.")
	parser.add_argument("-o", "--outfile", type=str, default="triallelic_by_sample.tab", help="outfile name")

	args = parser.parse_args()

	# vcf from stdin
	if sys.stdin.isatty():
		sys.exit("Error: Needs vcf from stdin, usage: 'zcat chr21.vcf.gz | count_triallelic_by_sample.py sample_file'")
	vcf = sys.stdin
	

	########

	# skip header until sample line
	vcf_header = D.get_sample_header(vcf)

	# get list of samples
	samples = D.get_pops_from_file(args.samples)
	
	# get col-numbers for every input sample
	sample_index = D.get_sample_colnumber(vcf_header, samples)
	
	# all other colnums with genotypes
	outside_index = [i for i in range(9, len(vcf_header)) if i not in sample_index.values()]

	# Count triallelic / biallelic sites
	sin, tri, tri_outside, bi_outside = count_triallelic_by_sample(vcf, sample_index, outside_index)
	
	# write outfile
	with open(args.outfile, "w") as out:
		for k,v in sin.items():
			out.write(k + '\t' + str(v) + '\n')
		out.write('trial_samples' + '\t' + str(tri) + '\n')
		out.write('trial_outside_samples' + '\t' + str(tri_outside) + '\n')
		out.write('bial_outside_samples' + '\t' + str(bi_outside) + '\n')



# CAUTION this assumes pure genotype-input
def count_triallelic_by_sample(vcf, sample_index, outside_index):
	'''
	sample_index: dict{sample:column-index}
	outside_index: list [column-indices] for all non-input samples  
	'''
	# initialize output
	print(sample_index)
	singletons = {}
	for s in sample_index.keys():
		singletons[s] = 0

	trial_samples = 0
	trial_outside_samples = 0
	bial_outside_samples = 0
	
	moreal = 0
	
	for line in vcf:
		try:
			line = line.rstrip().split()
			## bases/genotypes
			ref = line[3]
			alt = line[4]
			## remove invariable
			if alt == ".":
				continue
			## split ALT
			alt_alleles = alt.split(",")
			## skip more then 3 alleles
			if (len(alt_alleles) > 2):
				moreal += 1
				continue
			# remaining site are biallelic or triallelic
			
			## unique variants outside samples 
			outside_alleles = D.unique_variants(line, outside_index)
			
			## if biallelic: check if its biallelic outside samples
			if len(alt_alleles) == 1 and len(outside_alleles) == 2:
				bial_outside_samples += 1
			## if triallelic, check cases
			elif len(alt_alleles) == 2:
				# triallelic outside sample
				if len(outside_alleles) == 3:
					trial_outside_samples += 1
				# biallelic outside sample (dont count cases where one-all -> tri-all only by samples)
				elif len(outside_alleles) == 2:
					bial_outside_samples += 1
					# check if samples cause third state
					# (for perfect VCFs this would always be true, but here the ALT-field
					#  from SGDP is not (yet) updated after genotype-filtering, so that there
					#  are alleles in ALT which are not in any genotype)
					sample_alleles = D.unique_variants(line, list(sample_index.values()))
					if len(set(sample_alleles + outside_alleles)) == 3:
						trial_samples += 1
					## check if caused by singleton in one of the input-samples 
					# get sample alleles
					single_alleles = {}
					for k,v in sample_index.items():
						single_alleles[k] = D.unique_variants(line, [v])
					# check which ones are different from outside alleles
					causer = []
					for s,a in single_alleles.items():
						if len(set(a + outside_alleles)) > 2:
							causer.append(s)
					# (causer might be empty, ALT can have more alleles than what is really there in GTs)
					# check if its singleton
					if (len(causer) == 1):
						singletons[causer[0]] += 1
		except(IndexError, ValueError) as errore:
			sys.stderr.write(str(errore) + "\n")
			sys.stderr.write("\t".join(line) + "\n")
			sys.exit()
	sys.stdout.write('skipped '+ str(moreal) + ' sites with more than 3 alleles.\n')

	return singletons, trial_samples, trial_outside_samples, bial_outside_samples
	


if __name__ == "__main__":
	main()



