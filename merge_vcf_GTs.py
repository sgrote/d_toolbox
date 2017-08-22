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



def merge_alt_alleles(alt1, alt2, gt1_list, gt2_list):
	'''
	merge genotypes of two vcf-files when ALTs are not the same
	input: ALT-alleles-string from both files, GT-string lists from both files
	output: tuple of 3: new ALT, GT-string-lists from both files with new genotypes
	CAUTION: if gt1/gt2 have more than genotype, e.g. "1/0:23:xy", first entry is assumed to be the GT
	'''
	# check if one of the ALTs is missing and replace with the other one
	# no GT-altering needed
	if alt1 == ".":
		 return alt2, gt1_list, gt2_list
	if alt2 == ".":
		 return alt1, gt1_list, gt2_list
	# convert ALTs to list and merge 
	alt1 = alt1.split(",")
	alt2 = alt2.split(",")
	alt_merged = list(set(alt1 + alt2))
	alt_merged.sort()	
	# new ALT-numbers
	alt1_indi = [alt_merged.index(old_alt) + 1 for old_alt in alt1]
	alt2_indi = [alt_merged.index(old_alt) + 1 for old_alt in alt2]
	# new ALT-string
	alt_out = ",".join(alt_merged)
	# replace genotypes with new numbers
	gt1_new = [replace(gt, alt1_indi) for gt in gt1_list]
	gt2_new = [replace(gt, alt2_indi) for gt in gt2_list]
	return alt_out, gt1_new, gt2_new

''' test
merge_alt_alleles('A,C', 'C,T', ['0/1','2/1'], ['2|0','1|1']) == ('A,C,T', ['0/1', '2/1'], ['3|0', '2|2'])
merge_alt_alleles('A,C', 'C,A', ['0|1','1|2'], ['2/0','0/1']) == ('A,C', ['0|1', '1|2'], ['1/0', '0/2'])
merge_alt_alleles('C,A', 'A,C', ['0/1','1/2'], ['2/0','0/1']) == ('A,C', ['0/2', '2/1'], ['2/0', '0/1'])
merge_alt_alleles('A,C', 'AC,T', ['0/1','1/2'], ['1/0','0/2']) == ('A,AC,C,T', ['0/1', '1/3'], ['2/0', '0/4'])
merge_alt_alleles('A,C', 'AC,C', ['0/1','1/2'], ['1/0','0/2']) == ('A,AC,C', ['0/1', '1/3'], ['2/0', '0/3'])
merge_alt_alleles('G', 'T,G', ['0|1','0|0'], ['1|2','0|0']) == ('G,T', ['0|1', '0|0'], ['2|1', '0|0'])
merge_alt_alleles('G', 'T,G', ['0/1','0/0'], ['1/0','0/0']) == ('G,T', ['0/1', '0/0'], ['2/0', '0/0'])
merge_alt_alleles('T', '-,T', ['0/1','0/0'], ['1/0','0/2']) == ('-,T', ['0/2', '0/0'], ['1/0', '0/2'])
merge_alt_alleles('T', '-', ['0/1','0/0'], ['1/0','0/1']) == ('-,T', ['0/2', '0/0'], ['1/0', '0/1'])
merge_alt_alleles('T', 'T', ['0/1:15','0/0:12:13'], ['1/0']) == ('T', ['0/1:15', '0/0:12:13'], ['1/0'])
merge_alt_alleles('T', 'C', ['0/1','0/0'], ['1/0','0/1']) == ('C,T', ['0/2', '0/0'], ['1/0', '0/1'])
merge_alt_alleles('G', 'C', ['0/.','0/1'], ['.|.','./1']) == ('C,G', ['0/.', '0/2'], ['.|.', './1'])
merge_alt_alleles('T', 'C,T', ['1','0'], ['2','.']) == ('C,T', ['2', '0'], ['2', '.'])
merge_alt_alleles('AC,T,AAA', 'C,A,AC,-,G,GT,T,GTC,CCC,GA,TG', ['0/1','3/0'], ['10/0','7/1','11/4','9/3']) == ('-,A,AAA,AC,C,CCC,G,GA,GT,GTC,T,TG', ['0/4', '3/0'], ['8/0', '11/5', '12/1', '6/4'])
merge_alt_alleles('T', 'C,T', ['1:34:13','0:5'], ['2','.']) == ('C,T', ['2:34:13', '0:5'], ['2', '.'])
merge_alt_alleles('A', '.', ['0/1','1/1'], ['0|0','0|0']) == ('A', ['0/1', '1/1'], ['0|0', '0|0'])
merge_alt_alleles('.', 'A', ['0/0','0/0'], ['1|0','0|1']) == ('A', ['0/0', '0/0'], ['1|0', '0|1'])
merge_alt_alleles('.', '.', ['0/0','0/0'], ['0|0','0|0']) == ('.', ['0/0', '0/0'], ['0|0', '0|0'])

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
		alt, gt1, gt2 = merge_alt_alleles(line1[4], line2[4], line1[9:], line2[9:])
		out = line1[:4] + [alt] + line1[5:9] + gt1 + gt2
	out[5:8] = [".",".","."]
	return out

''' test
line1 = ['21','148','.','C','T','10','LowQual','AC=27','GT:GF','0/1:3','./0:-1']
line2 = ['21','148','.','C','A','10','PASS','.','GT','0|1','0|1']
line3 = ['21','148','.','C','.','10','PASS','.','GT','0|0','0|0']
line4 = ['21','148','.','C','A,T','10','PASS','.','GT','0|1','0|2']
line5 = ['21','148','.','C','G,A','10','PASS','.','GT','0|1','0|2']

merge_lines(line1, line2) == ['21','148','.','C','A,T','.','.','.','GT:GF','0/2:3','./0:-1','0|1','0|1']
merge_lines(line2, line1) == ['21','148','.','C','A,T','.','.','.','GT','0|1','0|1','0/2:3','./0:-1']
merge_lines(line3, line3) == ['21','148','.','C','.','.','.','.','GT','0|0','0|0','0|0','0|0']
merge_lines(line2, line2) == ['21','148','.','C','A','.','.','.','GT','0|1','0|1','0|1','0|1']
merge_lines(line1, line3) == ['21','148','.','C','T','.','.','.','GT:GF','0/1:3','./0:-1','0|0','0|0']
merge_lines(line1, line4) == ['21','148','.','C','A,T','.','.','.','GT:GF','0/2:3','./0:-1','0|1','0|2']
merge_lines(line4, line5) == ['21','148','.','C','A,G,T','.','.','.','GT','0|1','0|3','0|2','0|1']

'''




	

