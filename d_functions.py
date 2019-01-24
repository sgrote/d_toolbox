#!/usr/bin/python3

'''
Create all pairwise population matches
Use Info table to get population and sex for each individual
Go through vcf-header and collect all column-numbers that belong to distinct populations
Use this info to calculate abba, baba counts per block and population match

'''

#import StringIO
#import csv
import pandas as pd
import sys
import os.path


def get_range(chrom ,ranges_file):
	'''
	in: chrom (e.g. 21 or chr21),
		bed-file (chr, start, end) also 21 or chr21 (centromeres or dimensions)
	out: [start, end] of entry for chromosome
	'''
	# check for chr-prefix
	chrom = str(chrom)
	if not chrom.startswith("chr"):
		chrom = "chr" + chrom
	with open(ranges_file,"r") as ranges:
		for line in ranges:
			line = line.rstrip().split()
			line_chr = line[0]
			if not line_chr.startswith("chr"):
				line_chr = "chr" + line_chr
			if line_chr == chrom:
				return [int(line[1]), int(line[2])]
	print("Error, Chromosome range not found")


def get_pops_from_args(apop1, apop2, apop3, apop4):
	'''
	in: pop-arguments (either file-names or comma-separated strings)
	out: list [[pop1][pop2][pop3][pop4]] for single populations per 'column'
	'''
	if os.path.isfile(apop1): 
		print("Looks like pop-files.")
		return(get_pops_from_files(apop1, apop2, apop3, apop4))
	else:
		print("Looks like pop-lists.")
		poplists = [apop1, apop2, apop3, apop4]
		pops = list()
		for p in poplists:
			pops.append(p.split(","))
		return pops


def get_pops_from_files(pop1, pop2, pop3, pop4):
	'''
	in: files with one column containing population names
	out: list [[pop1][pop2][pop3][pop4]] for single populations per 'column'
	'''
	pop_files = [pop1, pop2, pop3, pop4]
	pops = [[] for k in range(4)]
	# read populations
	for i in range(4):
		with open(pop_files[i], "r") as p:
			for pop in p:
				pops[i].append(pop.rstrip())
	return pops



def get_pops_from_file(pop_file):
	'''
	in: file with one column containing population names
	out: list [pop1, pop2, ...]
	'''
	pops = []
	with open(pop_file, "r") as p:
		for pop in p:
			pops.append(pop.rstrip())
	return pops



def check_pw_pops(pops):
	'''
	check that pops [[pop1][pop2][pop3][pop4]] have same length
	'''
	l = len(pops[0])
	if any([len(pops[i]) != l for i in range(1,4)]):
		sys.exit("pop-lists have differing lengths (not valid for option -f/--fixedpairs). ")
	return True

''' test
check_pw_pops([["a","b"], ["c","d"], ["e","f"], ["s","f"]])
check_pw_pops([["a","b"], ["c","d"], ["e","f"], ["s"]])
'''

def get_pw_pops(pops):
	'''
	in: list [[pop1][pop2][pop3][pop4]] for single populations per 'column'
	out: list [[pop1][pop2][pop3][pop4]] for all matches,
	 avoiding self-matches and duplicates (including D(w,x,y,z), D(x,w,y,z))
	'''
	pw_pops = [[] for k in range(4)]
	# generate all possible combinations (That's ugly!)
	for i in range(len(pops[0])):
		for j in range(len(pops[1])):
			for k in range(len(pops[2])):
				for l in range(len(pops[3])):
					pw_pops[0].append(pops[0][i])
					pw_pops[1].append(pops[1][j])
					pw_pops[2].append(pops[2][k])
					pw_pops[3].append(pops[3][l])
	# remove unnecessary matches
	# find indices to remove
	remove = []
	pairs = []
	combs = [""]*4
	for i in range(len(pw_pops[0])):
		# avoid e.g. (Vidija, Vindija, Eskimo, Mbuti)
		if len(set([pw_pops[0][i],pw_pops[1][i],pw_pops[2][i],pw_pops[3][i]])) < 4:
			remove.append(i)
		# avoid e.g. (Altai,Vindija,...) and (Vindija,Altai,...), (...,chimp,orang) and (...,orang,chimp)
		else:
			combs[0] = pw_pops[0][i] + pw_pops[1][i] + pw_pops[2][i]+ pw_pops[3][i]
			combs[1] = pw_pops[1][i] + pw_pops[0][i] + pw_pops[2][i]+ pw_pops[3][i]
			combs[2] = pw_pops[0][i] + pw_pops[1][i] + pw_pops[3][i]+ pw_pops[2][i]
			combs[3] = pw_pops[1][i] + pw_pops[0][i] + pw_pops[3][i]+ pw_pops[2][i]
			if any([comb in pairs for comb in combs]):
				remove.append(i)
			else:
				pairs.append(combs[0])
	# remove
	for i in range(len(pw_pops)):
		for j in sorted(remove, reverse=True):
			del pw_pops[i][j]
	if len(pw_pops[1]) == 0:
		sys.exit("\nError: could not create pairwise population matches - check the pop input files.")
	return pw_pops


def get_pw_pops_from_file(pw_pops_file):
	'''
	in: csv-file with (at least) 4 columns named: pop1, pop2, pop3, pop4 defining the population matches
	out: list [[pop1][pop2][pop3][pop4]] for all matches,
	'''
	pw_pops_table = pd.read_csv(pw_pops_file)
	pop1 = list(pw_pops_table['pop1'])
	pop2 = list(pw_pops_table['pop2'])
	pop3 = list(pw_pops_table['pop3'])
	pop4 = list(pw_pops_table['pop4'])
	pw_pops = [pop1, pop2, pop3, pop4]
	return pw_pops

''' test
get_pw_pops_from_file("/mnt/expressions/steffi/D/infofiles/capture_eval_pairings_BY.csv")
'''

def get_unique_pops(pw_pops):
	'''
	in: list of lists for pairwise population matches
	out: list of unique population names
	'''
	unique_pops = set()
	for d_position in pw_pops:
		for pop in d_position:
			unique_pops.add(pop)
	return unique_pops

''' test
pw_pops = [['popA', 'popB'], ['popC', 'popA'], ['popC', 'popD'], ['chimp', 'chimp']]
get_unique_pops(pw_pops)
'''


def get_sample_header(vcf):
	'''
	in: open vcf-file
	out: list [last header line (with sample-names)]
	 as a side-effect the vcf's header is removed
	'''
	v = vcf.readline()
	while v[:6] != "#CHROM":
		v = vcf.readline()
	return v.rstrip().split()


# allow multiple pops per sample
def get_pop_colnumbers(info_file, vcf_header, pops):
	'''
	in: info_file csv-file with columns [sample, pop, sex]
		vcf_header as list
		pops as set of pops needed for D-stats
	out: dictionary {population:[colnumbers]}
		 dictionary {colnumber:sex}
	'''
	pop_colnums = {}
	col_gender = {}
	info = pd.read_csv(info_file)
	# for each donor-id from vcf-header, get population and sex if present in info-file
	for i in range(9, len(vcf_header)):
		if any(info['sample'] == vcf_header[i]):
			sample_info = info[info['sample'] == vcf_header[i]]
			# add col-number to population
			for pop in sample_info['pop']:
				# NEW: add only pops that are needed for D-stats
				if pop in pops:
					if pop not in pop_colnums:
						pop_colnums[pop] = [i]
					else:
						pop_colnums[pop].append(i)
			# add sex to col-number 
			# (first entry of sample_info['sex'], since this column is all the same for one sample)
			col_gender[i] = sample_info['sex'].iloc[0]
	return pop_colnums, col_gender

''' test
vcf_header = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','AltaiNeandertal','Vindija33.19','Denisova']
info_file = '/mnt/expressions/steffi/D/infofiles/example.csv'
pops = {'Vindija33.19','Denisova','Neandertal'}     # no Altai
popcol, colsex = get_pop_colnumbers(info_file, vcf_header, pops)
popcol == {'Neandertal': [9, 10], 'Denisova': [11], 'Vindija33.19': [10]}
colsex == {9: 'female', 10: 'female', 11: 'female'}
'''

# just the colnumber of a sample (useful for preprocessed allele-freq input)
def get_sample_colnumber(vcf_header, pops):
	'''
	in: vcf_header as list
		pops as set of pops needed for D-stats
	out: dictionary {population:colnumber}
	'''
	sample_colnum = {}
	for p in pops:
		sample_colnum[p] = vcf_header.index(p)
	return sample_colnum

''' test
vcf_header = ['first_fields','Altai','Vindija','Denisova', 'Altai']
get_sample_colnumber(vcf_header, ["Denisova", "Altai"]) == {'Denisova': 3, 'Altai': 1}
'''


# useful if non-input samples should get masked (after modifying REF/ALT for input)
def unique_pop_colnums(pop_colnums, n_samples):
	'''
	in: dict{pop:[cols]}, nr of samples in vcf
	out: [unique pop columns] , [non-input columns]
	'''
	pop_cols = set()
	for cols in pop_colnums.values():
		pop_cols.update(set(cols))
	pop_cols = sorted(list(pop_cols))
	non_pop_cols = [i for i in range(9, n_samples+9) if i not in pop_cols]
	return pop_cols, non_pop_cols

''' test
pop_colnums = {"a":[10,12,9,14], "b":[9,10,13]}
a, b = unique_pop_colnums(pop_colnums, 7)
a == [9,10,12,13,14]
b == [11,15]

'''


def unique_variants(vcf_line, vcf_colnums):
	'''
	return a set of uniqe allele indices
	from vcf_line[vcf_colnums]
	assumes genotypes only
	''' 
	gts = [vcf_line[i] for i in vcf_colnums]
	gt_string = "".join(gts)
	gt_chars = set(gt_string)
	gt_chars_clean = sorted([a for a in gt_chars if a not in "./|"])
	return gt_chars_clean 

''' test
vcf = ['first_fields', '0/3', '.|.', './1', '2|2', '0|.']
unique_variants(vcf, [1,5]) == ['0', '3']
unique_variants(vcf, [2,3,4]) == ['1', '2']
'''


def change_alt_index(vcf_line, vcf_colnums, old_al_index, new_al_index):
	'''
	modify input in place
	in all vcf_line[vcf_colnums] replace old_al_index with new_all_index
	either some-ALT -> REF, or some-ALT -> only ALT
	update REF or ALT as needed
	CAUTION: this requires that target index is NOT present in data
	'''
	# if old and new are the same, only update ALT if needed
	if old_al_index == new_al_index:
		if new_al_index == "0":
			return None
		if new_al_index == "1":
			vcf_line[4] = vcf_line[4].split(",")[0]
			return None
	if old_al_index == "0":
		sys.stderr.write("Why would a REF-allele be replaced?")
		sys.exit()
	# replace in genotype
	for i in vcf_colnums:
		vcf_line[i] = vcf_line[i].replace(old_al_index, new_al_index)
	# update REF or ALT
	new_al = vcf_line[4].split(",")[int(old_al_index)-1]
	if new_al_index == "1":
		vcf_line[4] = new_al # ALT
	elif new_al_index == "0":
		vcf_line[3] = new_al # REF
		# Dont update ALT here, to keep order for a second replacement!
		# (ALT -> REF will often have a following step, e.g. 1/3 -> 0/3 -> 0/1
	else:
		sys.stderr.write("New index must be 0 or 1: REF or one and only ALT")
		sys.exit()
	return None

''' test
line = ['21','148','.','C','A,T,G','0','.','.','GT','0/2','1/1','./0','2/.']
# nothing changes, C stays REF
change_alt_index(line, [9,10,11,12], "0", "0")
line == ['21','148','.','C','A,T,G','0','.','.','GT','0/2','1/1','./0','2/.']
# only update ALT, A stays first ALT
change_alt_index(line, [9,10,11,12], "1", "1")
line == ['21','148','.','C','A','0','.','.','GT','0/2','1/1','./0','2/.']
# T becomes only ALT in cols 9+11
line = ['21','148','.','C','A,T,G','0','.','.','GT','0/2','1/1','./0','2/.']
change_alt_index(line, [9, 11], "2", "1")
line == ['21', '148', '.', 'C', 'T', '0', '.', '.', 'GT', '0/1', '1/1', './0', '2/.']
# A becomes REF
line = ['21','148','.','C','A,T,G','0','.','.','GT','1|.','1|3']
change_alt_index(line, [9, 10], "1", "0")
line == ['21','148','.','A','A,T,G','0','.','.','GT','0|.','0|3']
# G becomes REF
line = ['21','148','.','C','A,T,G','0','.','.','GT','1|.','1|3']
change_alt_index(line, [9, 10], "3", "0")
line == ['21','148','.','G','A,T,G','0','.','.','GT','1|.','1|0']
# ALT is insertion
line = ['21','148','.','C','ATG','0','.','.','GT','1|.','1|0']
change_alt_index(line, [9, 10], "0", "0")
line == ['21','148','.','C','ATG','0','.','.','GT','1|.','1|0']
change_alt_index(line, [9, 10], "1", "1")
line == ['21','148','.','C','ATG','0','.','.','GT','1|.','1|0']

'''



def check_bial(vcf_line, pop_cols, non_pop_cols):
	'''
	check if a vcf-line is biallelic across	samples from pop_colnums ({pop:[col,nums,...]})
	return a modified line with one ALT-allele and {0,1} GT in input samples
	return None if its not bialleleic across samples
	(if line gets changed all non-input-samples become masked './.',
	 this not necessary here but would be very ugly if they stay unchanged)
	'''
	alt = vcf_line[4]
	# just to make sure
	if alt == ".":
		return None
	# unique allele-indices for input pops
	unique_allele_indis = unique_variants(vcf_line, pop_cols)
	if len(unique_allele_indis) != 2:
		return None
	if len(alt) == 1:
		return vcf_line # not check earlier because it could still be invariable across samples
	
	## update allele indices for input samples
	out = vcf_line[:]
	# replace REF first, then make the other ALT allele the one and only one
	# e.g. 1/3 -> 0/3 -> 0/1;  (unique alleles are sorted)
	change_alt_index(out, pop_cols, unique_allele_indis[0], "0")
	change_alt_index(out, pop_cols, unique_allele_indis[1], "1")
	# mask all other samples
	for i in non_pop_cols:
		out[i] = "./."
	
	return out

''' test
line = ['21','148','.','C','A,T,G','0','.','.','GT','0/2','1/1','./0','./.','1|.','1|3']
p1, np1 = [9, 10], [11,12,13,14] # trial
p2, np2  = [10,11,12,13], [9,14] # bial 1st allele
p3, np3 = [9,11], [10,12,13,14] # bial 2nd allele
p4, np4 = [10,12,13,14], [9,11] # bial 1st and 3rd alleles (modify REF!)

check_bial(line, p1, np1) == None
check_bial(line, p2, np2) == ['21','148','.','C','A','0','.','.','GT','./.','1/1','./0','./.','1|.','./.']
check_bial(line, p3, np3) == ['21','148','.','C','T','0','.','.','GT','0/1','./.','./0','./.','./.','./.']
check_bial(line, p4, np4) == ['21','148','.','A','G','0','.','.','GT','./.','0/0','./.','./.','0|.','0|1']

'''



def get_p(vcf_line, pop, digits=None, counts=False):
	'''
	alternative allele-freqs for one population for one biallelic site
	in: vcf = vcf-line as list
		pop = [colnumbers] for the population
		digits = digits to round to
		counts = return ("ALT/total") counts as string instead of the fraction 
	out: ALT-allele-freq
	'''
	
	ac = 0
	sum_alleles = 0
	
	## pre-check for all-REF or all missing
	genotypes = [vcf_line[colnumber] for colnumber in pop]
	if not counts and all([g in ["0/0", "0|0", "0"] for g in genotypes]):
		return 0.0
	if all([g in ["./.", ".|.", "."] for g in genotypes]):
		return None
	
	# count alternative alleles and total number of alleles for one population  
	for i in range(len(pop)):
		gt = genotypes[i]
		gt = gt.replace('/', '|')
		alleles = gt.split("|")
		for a in alleles:
			if a == "1":
				ac += 1
				sum_alleles += 1
			elif a == "0":
				sum_alleles += 1
			elif a != ".":
				sys.stderr.write("unexpected genotype\n")
				sys.stderr.write(" ".join(["genotype:", gt, ", allele:", a]) + "\n")
				sys.exit()
	
	# compute frequency 
	if sum_alleles == 0:
		return None   # unknown sites may be fully genotyped invariable sites, here there is simply no GT 
	elif counts:
		daf = str(ac) + "/" + str(sum_alleles)
	else: 
		daf = ac*1.0 / sum_alleles
		if digits:
			daf = round(daf, digits)
	return daf

''' test
vcf_line1 = ['21','148','.','C','T','0','.','.','GT','0/1','1/1','./0','./.','1|.','1|1','.|.','0/0','0|0']
vcf_line2 = ['X','148','.','C','T','0','.','.','GT','1/1','1','./0','./.','1|.','1/0','.','0/.','0|0']

get_p(vcf_line1, [9,10]) == 0.75
get_p(vcf_line1, [9,10,11]) == 3.0/5
get_p(vcf_line1, [9,10,11], counts=True) == '3/5'
get_p(vcf_line1, [9,10,11,14,15]) == 5.0/7
get_p(vcf_line1, [9,10,11,14,15], counts=True) == '5/7'
get_p(vcf_line1, [9,10,11,14,15], digits=4) == 0.7143
get_p(vcf_line1, [12,13,14,15]) == 1
get_p(vcf_line2, [9,10]) == 1
get_p(vcf_line2, [9,11]) == 2.0/3
get_p(vcf_line2, [9,14]) == 0.75
get_p(vcf_line2, [9,14], digits=1) == 0.8
get_p(vcf_line2, [10,14]) == 2.0/3
get_p(vcf_line2, [9,10,11,12]) == 0.75
get_p(vcf_line2, [9,10,11,12,13]) == 0.8
get_p(vcf_line1, [12,15]) == None
get_p(vcf_line2, [12,15]) == None
get_p(vcf_line2, [12,15], counts=True) == None
get_p(vcf_line2, [16,17]) == 0
get_p(vcf_line2, [16,17], counts=True) == '0/3'
get_p(vcf_line1, [16,17]) == 0
get_p(vcf_line1, [16,17], counts=True) == '0/4'

'''

	

# one list for containig block-abba-baba-sums for altai and vindija together in one list
# 0. blocks, 1. abba, baba, (bbba, aaba) 2. pw_pops  [ [[abba][baba]] [[abba][baba]] ...]
# CUATION this assumes pure genotype-input
def abba_block_sums(vcf, pop_colnums, pw_pops, block_size, centro_range=None, transver=False, sitesfile=None, af_input=False, min_p3=0, max_p3=1, count_bbba=False):
	pops = set(pw_pops[0]+pw_pops[1]+pw_pops[2]+pw_pops[3])
	if sitesfile:
		sites_file = open("sites","w")
		if sitesfile == "sites":
			sites_file.write('\t'.join(["block","pos"]) + "\n")
		elif sitesfile == "full":
			sites_file.write('\t'.join(["block","pos","ref","alt"] + list(pops)) + "\n")
	# initialize sums
	block_count = 0
	#n_blocks = len(block_borders)-1
	blocksites = 0
	n_comp = len(pw_pops[0])
	# initialize abba and baba block-sums [[abba],[baba], (bbba, aaba)]
	if count_bbba:
		block_sum = [[0]*n_comp, [0]*n_comp, [0]*n_comp, [0]*n_comp]
	else:
		block_sum = [[0]*n_comp, [0]*n_comp]
	sites_comp = [0]*n_comp
	blocks = []
	trans_county = 0
	nbial_county = 0
	# convert blocksize from kb to bases
	block_size = block_size * 1000
	first = True
	for line in vcf:
		line = line.rstrip().split()
		# skip invariable sites
		if line[4] == ".":
			continue
		# get first pos for block-border
		if first:
			block_end = int(line[1]) + block_size
			print("starting block 1 until position", block_end)
			first = False
		# exclude centromere (usually not needed because of manifesto)
		if centro_range and (int(line[1]) > centro_range[0] and int(line[1]) < centro_range[1]):
			#print("Skipping centromere position", int(line[1]), centro_range)
			continue
		## bases/genotypes
		ref = line[3]
		alt = line[4]
		## reduce to biallelic 
		if alt not in ["A","C","T","G"] or ref not in ["A","C","T","G"]:
			#print("Skipping non-biallelic:", line[:7])
			nbial_county += 1
			continue
		## reduce to transversions 
		if transver: 
			if ("A" in ref+alt and "G" in ref+alt) or ("C" in ref+alt and "T" in ref+alt):
				#print("Skipping transition:", line[:7])
				trans_county += 1
				continue
		## next block?
		while int(line[1]) > block_end:
			# last entry, number of sites for the current block
			print(blocksites)
			# add current block and initialize again
			blocks.append(block_sum)
			if count_bbba:
				block_sum = [[0]*n_comp, [0]*n_comp, [0]*n_comp, [0]*n_comp]
			else:
				block_sum = [[0]*n_comp, [0]*n_comp]
			block_end = block_end + block_size
			blocksites = 0
			block_count += 1
			print("starting block ", block_count+1, "until position", block_end)
		
		## get alternative allele freqs for input populations
		p = {}
		## if allele-freq-input: get_p directly from file with pop_columns[pop]
		if af_input:
			for pop in pops:
				af = line[pop_colnums[pop]]
				p[pop] = None if af in ["None", "."] else float(af)
		else:
			for pop in pops:
				p[pop] = get_p(line, pop_colnums[pop])
		
		#print(p)
		## check that not all pop-freqs are 0 or 1 and None (pop with variant might not be part of D-stats)
		pset = set(p.values())
		if pset.issubset({0,None}) or pset.issubset({1,None}):
			#print("Uninformative site:", pset) 
			continue
		
		## compute p(ABBA) and p(BABA)  
		abba = []
		baba = []
		if count_bbba:
			bbba = []
			aaba = []
		for i in range(n_comp):
			p1 = p[pw_pops[0][i]]
			p2 = p[pw_pops[1][i]]
			p3 = p[pw_pops[2][i]]
			p4 = p[pw_pops[3][i]]
			# check if valid, else append 0
			if None in [p1, p2, p3, p4]:
				abba.append(0)
				baba.append(0)
				if count_bbba:
					bbba.append(0)
					aaba.append(0)
				continue
			if min_p3 != 0 or max_p3 != 1:
				# convert from ALT-freq to derived B-freq given outgroup
				if p4 == 0:
					p3_derived = p3
				elif p4 == 1:
					p3_derived = 1 - p3 
				else:
					print(line)
					print("p(",pw_pops[3][i],") = ", p4)
					sys.exit("outgroup pop4 has to be fixed to restrict on derived allele freq in pop3 (to infer the derived state).")
				if p3_derived < min_p3 or p3_derived > max_p3:
					abba.append(0)
					baba.append(0)
					if count_bbba:
						bbba.append(0)
						aaba.append(0)
					continue
			# if valid compute ABBA, BABA
			p_abba = (p1*(1-p2)*(1-p3)*p4) + ((1-p1)*p2*p3*(1-p4))
			p_baba = ((1-p1)*p2*(1-p3)*p4) + (p1*(1-p2)*p3*(1-p4))
			abba.append(p_abba)
			baba.append(p_baba)
			if count_bbba:
				p_bbba = ((1-p1)*(1-p2)*(1-p3)*p4) + (p1*p2*p3*(1-p4))
				p_aaba = (p1*p2*(1-p3)*p4) + ((1-p1)*(1-p2)*p3*(1-p4))
				bbba.append(p_bbba)
				aaba.append(p_aaba)
			# add to sites-count-per pw_pop if informative
			if p_abba > 0 or p_baba > 0:
				#print(i)
				#print(sites_comp)
				sites_comp[i] += 1
		
		# add bbba and aaba to block_sum
		if count_bbba and not (sum(bbba)==0 and sum(aaba)==0):
			for i in range(n_comp):
				block_sum[2][i] += bbba[i]
				block_sum[3][i] += aaba[i]
		
		# check that not all ABBA BABA is 0
		if sum(abba)==0 and sum(baba)==0:
			#print("uninformative site", line)
			continue
		# add abba and baba to block_sum
		for i in range(n_comp):
			block_sum[0][i] += abba[i]
			block_sum[1][i] += baba[i]  
		## add to sites counter and sites file
		blocksites += 1
		if sitesfile and sitesfile == "sites":
			sites_file.write('\t'.join([str(block_count+1)] + [line[1]]) + "\n")
		elif sitesfile and sitesfile == "full":
			ps = [p[pop] for pop in  pops]
			sites_file.write("\t".join([str(block_count+1)] + [line[1]] + line[3:5]) + "\t" + "\t".join(map(str,ps)) + "\n")
		## NEW: add position and ALT-allele-freqs to sites files
	## save last block  
	blocks.append(block_sum)
	if sitesfile:
		sites_file.close()
	print("Number of skipped non-biallelic sites: %d" % nbial_county)
	print("Number of skipped transitions: %d" % trans_county)
	return blocks, sites_comp
	

	
## output formatting    

# create "dataframe" [[pop1][pop2][pop3][pop4][block][abba][baba]]  
# input: pw_pops [[pop1][pop2][pop3][pop4]; blocks [ [[abba][baba]] [[abba][baba]] ...]
def rearrange_blocks(pw_pops, blocks):  
	out = [[] for i in range(7)]
	n_blocks = len(blocks)
	n_comps = len(pw_pops[0])
	for i in range(len(pw_pops)):
		out[i] += pw_pops[i] * n_blocks
	for i in range(n_blocks):
		out[4] += [i+1] * n_comps
		out[5] += blocks[i][0]  # abba
		out[6] += blocks[i][1]  # baba
	return out

# same with additional BBBA and AABA counts
def rearrange_blocks_bbba(pw_pops, blocks):  
	out = [[] for i in range(9)]
	n_blocks = len(blocks)
	n_comps = len(pw_pops[0])
	for i in range(len(pw_pops)):
		out[i] += pw_pops[i] * n_blocks
	for i in range(n_blocks):
		out[4] += [i+1] * n_comps
		out[5] += blocks[i][0]  # abba
		out[6] += blocks[i][1]  # baba
		out[7] += blocks[i][2]  # bbba
		out[8] += blocks[i][3]  # aaba
	return out

# print a list of lists with same lengths like a data frame         
def print_table_to_file(table,ofile, mode="w"):
	with open(ofile, mode) as out:
		for i in range(len(table[0])):
			line = []
			for j in range(len(table)):
				line.append(table[j][i]) 
			line.append("\n")   
			out.write("\t".join(map(str, line)))

# print a list of lists with same lengths like a data frame         
def print_table(table):
	for i in range(len(table[0])):
		line = []
		for j in range(len(table)):
			line.append(table[j][i]) 
		print("\t".join(map(str, line)))
