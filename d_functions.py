#!/usr/bin/python

'''
Use Info table to get population for each individual
Go through vcf-header and collect all column-numbers that belong to distinct populations
Create all pairwise population matches
Use this info to calculate d-statistics on headerless paste of Nean-Denis-Chimp-Pongo-bases and human-vcf-file

ouput: [[pop1][pop2][abba_nean][baba_nean][abba_denis][baba_denis][#vindija_sites][#altai_sites]
'''

import StringIO
import csv

# extract genotype
# (e.g. ["0/1:50", "0/0:10"] -> ["0|1","0|0"])
def extract_gt(gts):
	gt_line = []
	for gt in gts:
		gt = gt.replace("/","|")
		gt_line.append(gt[:3])
	return(gt_line)

# NEW: account for indivdual quality "0/1:50:.." -> "0|1", "0/1:-1" -> ".|."
def extract_c_gt(c_gts):
	gt_line = []
	for gt in c_gts:
		genotype = gt.split(":")[0]
		qual = gt.split(":")[1]
		if int(qual) >= 0:
			genotype = genotype.replace("/","|")
			gt_line.append(genotype)
		else:
			gt_line.append(".|.")
	return(gt_line)

# [start, end] of centromere or chromosome
# input: chrom, start, end
def get_range(chrom ,ranges_file):
	with open(ranges_file,"r") as ranges:
		for line in ranges:
			line = line.rstrip()
			line = line.split()
			if line[0] == chrom:
				return([int(line[1]), int(line[2])]) 
	print("Error, Chromosome range not found")

# [borders] for chromosome blocks
# input-chrom_start, chrom_end (range of potentially informative sites from nean_denis_derived)
# input-centro_range: chrom, start, end (from UCSC cyto bands)
# if kb=True n_blocks defines the size of blocks instead of the number of blocks
def get_block_borders(chrom_start, chrom_end, centro_start, centro_end, n_blocks, kb="false"):	
	skipped = False
	n_blocks = int(n_blocks)
	chrom_start = int(chrom_start)
	chrom_end = int(chrom_end)
	chromsize = chrom_end - chrom_start + 1
	centrosize = centro_end - centro_start + 1
	## get effective chrom size (interesting part) without centromer
	# centromer lies completely in region
	if chrom_start <= centro_start and chrom_end >= centro_end:
		eff_chromsize = chromsize - centrosize
	# centromer overlaps at beginning	
	elif chrom_start > centro_start and chrom_start < centro_end:
		eff_chromsize = chromsize - (centro_end - chrom_start) 
		chrom_start += centro_end - chrom_start
		print "overlap at beginning, new chromsize=", eff_chromsize, chrom_start, chrom_end
		skipped = True	
	# centromer overlaps at end
	elif chrom_end > centro_start and chrom_end < centro_end:
		eff_chromsize = chromsize - (chrom_end - centro_start) 	
		skipped = True
	# centromer lies outside
	elif chrom_start > centro_end or chrom_end < centro_start:
		eff_chromsize = chromsize
		skipped = True
	if kb == "true": 	
		blocksize = n_blocks * 1000
		n_blocks =  eff_chromsize / blocksize
	else: 
		blocksize = eff_chromsize / n_blocks		
	# initialize with start of chromosome
	borders = [chrom_start - 1]
	for i in range(n_blocks):
		next_border = borders[-1] + blocksize
		# skip centromere (add centromere-size to block, centro positions are skipped in d_block_sums())
		# also ensure that start of informative positions lies before end of centromere
		if skipped==False and borders[0] <= centro_end and next_border >= centro_start:
			print "landed in centromere", next_border
			next_border += centrosize
			skipped = True
		borders.append(next_border)
	# for kb-windows: append window rather then extending the last one	
	if kb == "true": borders.append(borders[-1] + blocksize)
	# last block might be longer, but less than n_blocks (int division) for 
	else: borders[-1] = borders[-1] + n_blocks	
	return(borders)	


# dictionary {individual:[superpop, population, sex]}
def get_pop_dict(info_file):
	with open(info_file, 'r') as f:
		reader = csv.reader(f)
		reader.next()
		pop_dict={}
		# get superpop, population and gender for every individual
		for line in reader:			
			pop_dict[line[0]] = [line[2], line[3], line[4]]		
	return(pop_dict)


# dictionary {population:[colnumbers]}
# dictionary {colnumber:gender}
def get_pop_colnumbers(pop_dict, vcf_header):
	pop_colnums={}
	col_gender={}			
	# for each donor-id from vcf-header, get population
	for i in range(len(vcf_header)):
		if vcf_header[i] in pop_dict:
			# add col-number to population
			pop = pop_dict[vcf_header[i]][1]
			if pop not in pop_colnums:
				pop_colnums[pop] = [i]
			else:
				pop_colnums[pop].append(i)
			# add gender to col-number
			gender = pop_dict[vcf_header[i]][2]
			col_gender[i] = gender 	
			# NEW: Africa2
			if (pop_dict[vcf_header[i]][0]=="Africa") and (not pop in ["Somali", "Sahrawi", "Mozabite", "Masai"]):
				if "Africa2" not in pop_colnums:
					pop_colnums["Africa2"] = [i]
				else:
					pop_colnums["Africa2"].append(i)
	return(pop_colnums, col_gender)


# population-matches [[pop1][pop2][pop3]] pop1=all,including archs; pop2=mbuti,bantu,esan; 3=altai,vindija 
def get_pops_from_files(pop1, pop2, pop3, pop4):
	pop_files = [pop1, pop2, pop3, pop4]
	pops = [[] for k in range(4)]
	pw_pops = [[] for k in range(4)]
	# read populations
	for i in range(4):
		with open(pop_files[i], "r") as p:
			for pop in p:
				pops[i].append(pop.rstrip()) 
	print pops
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
	#remove = [i for i in range(len(pw_pops[0])) if len(set([pw_pops[0][i],pw_pops[1][i],pw_pops[2][i]])) < 3]
	# find indices to remove
	remove = []
	pairs = []
	for i in range(len(pw_pops[0])):
		# avoid e.g. (Vidija, Vindija, Eskimo, Mbuti)
		if len(set([pw_pops[0][i],pw_pops[1][i],pw_pops[2][i]])) < 3:
			remove.append(i)
		# avoid e.g. (Altai, Vindija, ...) and (Vindija, Altai, ...)
		else:
			comb1 = pw_pops[0][i] + pw_pops[1][i] + pw_pops[2][i]+ pw_pops[3][i]
			comb2 = pw_pops[1][i] + pw_pops[0][i] + pw_pops[2][i]+ pw_pops[3][i]
			if comb1 in pairs or comb2 in pairs:
				remove.append(i)
			else:
				pairs.append(comb1)
	# remove
	for i in range(len(pw_pops)):
		for j in sorted(remove, reverse=True):
			del pw_pops[i][j]
	return(pw_pops)


# get p (derived allele frequency) for a population given vcf and col_nums and derived-GT
# (b is either REF or ALT)
# requires 1|0 GT instead of 1/0
def get_p(vcf, pop, col_gender, X):
	ac = 0
	sum_alleles = 0
	# count alternative alleles and total number of alleles for one population	
	for indi in pop:
		# X-chromosome for males (handle 1|., 0, ., .|., 1|1, and so on)
		if X==True and col_gender[indi]=="male":
			# alternative alleles
			if "1" in vcf[indi]: ac += 1
			# all alleles 
			if not(vcf[indi] == "." or vcf[indi] == ".|."):
				sum_alleles += 1
		else:
			# alternative alleles
			if vcf[indi] == "1|1": ac += 2
			elif "1" in vcf[indi]: ac += 1
			# all alleles
			if "." not in vcf[indi]: sum_alleles += 2
			elif vcf[indi] != ".|.": sum_alleles += 1
	# compute frequency 
	if sum_alleles == 0:
		#print "This Pop has no genotypes:", pop
		return(None)   # unknown sites may be fully genotyped invariable sites, here there is simply no GT 
	else: daf = ac*1.0 / sum_alleles		

	return(daf)

	

# one list for containig block-abba-baba-sums for altai and vindija together in one list
# 0. blocks, 1. abba, baba, 2. pw_pops  [ [[abba][baba]] [[abba][baba]] ...]
def abba_block_sums(vcf, pop_colnums, pw_pops, col_gender, block_borders, centro_range, private="false", transver="false", affixed="false"):
	print "affixed:", affixed
	pops = set(pw_pops[0]+pw_pops[1]+pw_pops[2]+pw_pops[3]+["Africa2"])
	sites_file = open("sites","w")
	sites_file.write('\t'.join(["block","pos","ref","alt"] + list(pops)) + "\n")		
	# initialize sums
	block_count = 0
	n_blocks = len(block_borders)-1
	sites = [0]*n_blocks
	n_comp = len(pw_pops[0])
	# initialize abba and baba block-sums [[abba],[baba]] 
	block_sum = [[0]*n_comp, [0]*n_comp]
	sites_comp = [0]*n_comp
	blocks = []			
	trans_county = 0			
	affixed_county = 0			
	for line in vcf:
		line = line.rstrip().split()
		#print line
		## bases/genotypes
		ref = line[3]
		alt = line[4]		
		## reduce to transversions 
		if transver == "true": 
			if ("A" in ref+alt and "G" in ref+alt) or ("C" in ref+alt and "T" in ref+alt):
				#print "Skipping transition:", line[:7]
				trans_county += 1
				continue
		## exclude centromere, NEW: not needed because of manifesto and everything
		#if int(line[1]) > centro_range[0] and int(line[1]) < centro_range[1]:
			#print "Skipping centromere position", int(line[1]), centro_range
			#continue
		## next block?
		while int(line[1]) > block_borders[block_count+1]: 
			# add current block and initialize again
			blocks.append(block_sum)
			block_sum = [[0]*n_comp, [0]*n_comp]			
			print sites[block_count]
			block_count += 1
			print "".join(map(str,["starting block ", block_count+1]))
		## is it X-chrom?
		if line[0]=="X": X = True
		else: X = False	
		# get alternative allele freqs for all populations (NEW: only input populations)	
		p = {}
		for pop in pops:	
			p[pop] = get_p(line, pop_colnums[pop], col_gender, X)	
		#print p
		
		## NEW: reduce to almost fixed in Africa2 (0.02 threshold = 1 allele in 37 Africa2-indis)
		if affixed == "true" and not (p["Africa2"]<0.02 or p["Africa2"]>0.98):
			#print p["Africa2"]
			affixed_county += 1
			continue
		## TODO: apply this somehow to new 'derived' situation?
		## skip if sites are required to be derived only in Neandertals but not in Denisovan
		#if private=="true" and p["Denisovan"] !=0:
			#continue
			
		# compute p(ABBA) and p(BABA) 	
		abba = []
		baba = []
		for i in range(n_comp):
			#if affixed == "true":
				#if (p[pw_pops[2][i]]==1 and p["Africa2"]>0.02) or (p[pw_pops[2][i]]==0 and p["Africa2"]<0.98):
					#abba.append(0)
					#baba.append(0)
					#continue
			if None in [p[pw_pops[0][i]], p[pw_pops[1][i]], p[pw_pops[2][i]], p[pw_pops[3][i]]]:
				abba.append(0)
				baba.append(0)
			else:
				p1 = p[pw_pops[0][i]]
				p2 = p[pw_pops[1][i]]
				p3 = p[pw_pops[2][i]]
				p4 = p[pw_pops[3][i]]
				p_abba = (p1*(1-p2)*(1-p3)*p4) + ((1-p1)*p2*p3*(1-p4))
				p_baba = ((1-p1)*p2*(1-p3)*p4) + (p1*(1-p2)*p3*(1-p4))
				abba.append(p_abba)
				baba.append(p_baba)
				# NEW add to sites-count-per pw_pop if informative
				if p_abba > 0 or p_baba > 0:
					#print i
					#print sites_comp				
					sites_comp[i] += 1
		# check that not all is 0
		if sum(abba)==0 and sum(baba)==0:
			#print "uninformative site", line
			continue 
		# add abba and baba to block_sum
		for i in range(n_comp):
			block_sum[0][i] += abba[i]
			block_sum[1][i] += baba[i]	
		# NEW: add position and ALT-allele-freqs to sites files
		ps = [p[pop] for pop in  pops]
		sites_file.write("\t".join([str(block_count+1)] + [line[1]] + line[3:5] + map(str,ps)) + "\n")
		sites[block_count] += 1
		if (sum(sites) % 200 == 0):
			print sum(sites)			
	# save last block	
	blocks.append(block_sum)
	sites_file.close()
	print "Number of skipped transitions: %d" % trans_county							
	print "Number of skipped sites too variable in Africans: %d" % affixed_county						
	return(blocks, sites_comp)
	

	
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
		out[5] += blocks[i][0] 	# abba
		out[6] += blocks[i][1] 	# baba
	return(out)


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
