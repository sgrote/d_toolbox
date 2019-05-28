#!/usr/bin/python3

import string

def get_new_indi(al, new_indi):
	'''
	convert allele to new alt-number if int and >0 ; if not return input ("0",".")
	'''
	if al=="0":
		return al
	try:
		return new_indi[int(al)-1]
	except ValueError:
		return al   


def replace(gt, new_indi):
	'''
	replace every number i in a GT with element i-1 in the new index for the old ALT alleles
	leave "." as "." and "0" as "0" and account for 1/0:23:45 etc = gt:other:info
	'''
	if("|" in gt): sep = "|"
	elif("/" in gt): sep = "/"
	else: sep = ""
	fields = gt.split(":")
	gt = fields[0]
	if not sep == "":
		gt = gt.split(sep)		
	new_alleles = [get_new_indi(al, new_indi) for al in gt]
	new_alleles = sep.join(map(str,new_alleles))
	if len(fields) > 1:
		appendix = fields[1:]
		new_alleles = ":".join([new_alleles] + appendix)
	return new_alleles

''' test
replace("1|2",[3,2,1]) == "3|2"
replace("1|2:0:12",[3,2,1]) == "3|2:0:12"
replace("./2",[3,2,1]) == "./2"
replace(".|.",[3,2,1]) == ".|."
replace(".",[3,2,1]) == "."
replace("1",[3,2,1]) == "3"
replace("3",[3,2,1]) == "1"
replace("2",[3,2,1]) == "2"
replace("0",[3,2,1]) == "0"
replace("1/2:12:13",[3,2,1]) == "3/2:12:13"
replace(".|.:0",[3,2,1]) == ".|.:0"
replace("3",[3,2,1]) == "1"

'''



def merge_alt_alleles(alt1, alt2, gt2_list):
	'''
	merge genotypes of two vcf-files when ALTs are not the same
	input: ALT-alleles-string from both files, GT-string list from second file
	output: tuple of 2: new ALT, GT-string-list for second file with new genotypes
	CAUTION: if gt1/gt2 have more than genotype, e.g. "1/0:23:xy", first entry is assumed to be the GT
	'''
	# check if one of the ALTs is missing and replace with the other one
	# no GT-altering needed
	if alt1 == ".":
		 return alt2, gt2_list
	if alt2 == ".":
		 return alt1, gt2_list
	# convert ALTs to list and add ALTs from alt2 that are missing in alt1 
	alt1 = alt1.split(",")
	alt2 = alt2.split(",")
	alt_merged = alt1 + [a for a in alt2 if a not in alt1]
	# new ALT-numbers
	alt2_indi = [alt_merged.index(old_alt) + 1 for old_alt in alt2]
	# new ALT-string
	alt_out = ",".join(alt_merged)
	# replace genotypes with new numbers
	gt2_new = [replace(gt, alt2_indi) for gt in gt2_list]
	return alt_out, gt2_new

''' test
merge_alt_alleles('A,C', 'C,T', ['2|0','1|1']) == ('A,C,T', ['3|0', '2|2'])
merge_alt_alleles('A,C', 'C,A', ['2/0','0/1']) == ('A,C', ['1/0', '0/2'])
merge_alt_alleles('C,A', 'A,C', ['2/0','0/1']) == ('C,A', ['1/0', '0/2'])
merge_alt_alleles('A,C', 'AC,T', ['1/0','0/2']) == ('A,C,AC,T', ['3/0', '0/4'])
merge_alt_alleles('A,C', 'AC,C', ['1/0','0/2']) == ('A,C,AC', ['3/0', '0/2'])
merge_alt_alleles('G', 'T,G', ['1|2','0|0']) == ('G,T', ['2|1', '0|0'])
merge_alt_alleles('G', 'T,G', ['1/0','0/0']) == ('G,T', ['2/0', '0/0'])
merge_alt_alleles('T', '-,T', ['1/0','0/2']) == ('T,-', ['2/0', '0/1'])
merge_alt_alleles('T', '-', ['1/0','0/1']) == ('T,-', ['2/0', '0/2'])
merge_alt_alleles('T', 'T', ['1/0']) == ('T', ['1/0'])
merge_alt_alleles('T', 'C', ['1/0','0/1']) == ('T,C', ['2/0', '0/2'])
merge_alt_alleles('G', 'C', ['.|.','./1']) == ('G,C', ['.|.', './2'])
merge_alt_alleles('T', 'C,T', ['2','.']) == ('T,C', ['1', '.'])
merge_alt_alleles('AC,T,AAA', 'C,A,AC,-,G,GT,T,GTC,CCC,GA,TG', ['10/0','7/1','11/4','9/3']) == ('AC,T,AAA,C,A,-,G,GT,GTC,CCC,GA,TG', ['11/0', '2/4', '12/6', '10/1'])
merge_alt_alleles('T', 'C,T', ['2','.']) == ('T,C', ['1', '.'])
merge_alt_alleles('A', '.', ['0|0','0|0']) == ('A', ['0|0', '0|0'])
merge_alt_alleles('.', 'A', ['1|0','0|1']) == ('A', ['1|0', '0|1'])
merge_alt_alleles('.', '.', ['0|0','0|0']) == ('.', ['0|0', '0|0'])

'''

def merge_lines(line1, line2):
	'''
	merge two vcf-file lines
	input: lists of vcf-line fields
	output: list with merged line; QUAL, FILTER, INFO are replaced by "."; other is taken from line1
	CAUTION: if gt1/gt2 have more than genotype, e.g. "1/0:23:xy", first entry is assumed to be the GT
	'''
	if line1[4] == line2[4]:
		out = line1 + line2[9:]
	else: 
		# merge Alt-alleles CAUTION: this assumes genotype to be first entry of samples
		alt, gt2 = merge_alt_alleles(line1[4], line2[4], line2[9:])
		out = line1[:4] + [alt] + line1[5:] + gt2
	out[5:8] = [".",".","."]
	return out

''' test
line1 = ['21','148','.','C','T','10','LowQual','AC=27','GT:GF','0/1:3','./0:-1']
line2 = ['21','148','.','C','A','10','PASS','.','GT','0|1','0|1']
line3 = ['21','148','.','C','.','10','PASS','.','GT','0|0','0|0']
line4 = ['21','148','.','C','A,T','10','PASS','.','GT','0|1','0|2']
line5 = ['21','148','.','C','G,A','10','PASS','.','GT','0|1','0|2']

merge_lines(line1, line2) == ['21','148','.','C','T,A','.','.','.','GT:GF','0/1:3','./0:-1','0|2','0|2']
merge_lines(line2, line1) == ['21','148','.','C','A,T','.','.','.','GT','0|1','0|1','0/2:3','./0:-1']
merge_lines(line3, line3) == ['21','148','.','C','.','.','.','.','GT','0|0','0|0','0|0','0|0']
merge_lines(line2, line2) == ['21','148','.','C','A','.','.','.','GT','0|1','0|1','0|1','0|1']
merge_lines(line1, line3) == ['21','148','.','C','T','.','.','.','GT:GF','0/1:3','./0:-1','0|0','0|0']
merge_lines(line1, line4) == ['21','148','.','C','T,A','.','.','.','GT:GF','0/1:3','./0:-1','0|2','0|1']
merge_lines(line4, line5) == ['21','148','.','C','A,T,G','.','.','.','GT','0|1','0|2','0|3','0|1']
merge_lines(line5, line4) == ['21','148','.','C','G,A,T','.','.','.','GT','0|1','0|2','0|2','0|3']

'''

# this can be useful when individual genotypes were filtered, uses functions from here
def check_complete_alt(line, pop_cols=[], non_pop_cols=[]):
	'''
	check if all ALT-bases are really in data (optionally restricted for pop_cols)
	in: vcf_line as list of strings, optionally lists of target and non-target samples
	out: original or updated vcf_line
	(if line gets changed all non-input-samples become masked './.',
	 this not necessary here but would be very ugly if they stay unchanged)
	'''
	alt = line[4]
	if alt == ".":
		return line
	alt_alleles = alt.split(",")
	# restrict to specific samples?
	if len(pop_cols) == 0:
		pop_cols = [i for i in range(9, len(line))]
	# extract genotype from more complex field (first entry)
	gts = [line[i] for i in pop_cols]
	gts_string = "".join([g.split(":")[0] for g in gts])
	presente = [str(i+1) in gts_string for i in range(len(alt_alleles))]
	# return input if all alleles are there
	if all(presente):
		return line
	# create new ALT field
	new_alt = [a for (a,p) in zip(alt_alleles, presente) if p]
	# replace by "." if none of the ALT alleles is really there
	if len(new_alt) == 0:
		alt_out = "."
		# no GT-modification needed, must be {0,.}
		out_line = line[:4] + [alt_out] + line[5:]
	else:
		alt_out = ",".join(new_alt)
		# update GTs (index of old alleles in new alleles order, None for removed alleles)
		new_alt_indi = [new_alt.index(o) + 1 if o in new_alt else None for o in alt_alleles]
		out_line = line
		out_line[4] = alt_out
		for i in range(len(pop_cols)):
			out_line[pop_cols[i]] = replace(gts[i], new_alt_indi)
		# mask other samples if any
		for i in non_pop_cols:
			out_line[i] = "./."
		
	return out_line


''' test
line0 = ['21','148','.','C','T,G','10','LowQual','AC=27','GT:GF','0/1:3','./2:-1'] # all ALT there
line1 = ['21','148','.','C','T,G','10','LowQual','AC=27','GT:GF','0/1:3','./0:-1'] # second ALT miss
line2 = ['21','148','.','C','T,G','10','LowQual','AC=27','GT:GF','0/2:3','./0:-1'] # first ALT miss
line3 = ['21','148','.','C','.','10','PASS','.','GT','0|0','0|0'] # no ALT at all
line4 = ['21','148','.','C','G,A','10','PASS','.','GT','0|0','0|0'] # all ALT miss
line5 = ['21','148','.','C','G,A,T','10','PASS','.','GT','0|1','2|0','0|3','1|2'] # more samples

check_complete_alt(line0) == ['21','148','.','C','T,G','10','LowQual','AC=27','GT:GF','0/1:3','./2:-1']
check_complete_alt(line0, [9], [10]) == ['21','148','.','C','T','10','LowQual','AC=27','GT:GF','0/1:3','./.']
check_complete_alt(line1) == ['21','148','.','C','T','10','LowQual','AC=27','GT:GF','0/1:3','./0:-1']
check_complete_alt(line2) == ['21','148','.','C','G','10','LowQual','AC=27','GT:GF','0/1:3','./0:-1']
check_complete_alt(line3) == ['21','148','.','C','.','10','PASS','.','GT','0|0','0|0']
check_complete_alt(line4) == ['21','148','.','C','.','10','PASS','.','GT','0|0','0|0']
check_complete_alt(line5, [9,11], [10,12]) == ['21','148','.','C','G,T','10','PASS','.','GT','0|1','./.','0|2','./.']

'''
