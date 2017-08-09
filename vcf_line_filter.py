
'''
filter_with_string:
take a vcf-line and a ONE string with a filter expression per line
e.g. "(QUAL > 0 or FILTER != 'LowQual') and AC > 20"
may contain keywords QUAL, FILTER and all keywords from info (to allow any combination of and/or)
return True/False
CAUTION: string with filter-expression needs whitespaces, at least around keywords and operators
(this is ugly but independent from keyword-info in vcf-header, which might be missing or incomplete)
'''

'''
do filtering for individual genotypes in a separate process
'''


def get_keys(filter_string):
	'''
	in: list with condition expression
	out: keywords and their positions (before relational operator)
	'''
	key_indis = [i-1 for i, x in enumerate(filter_string) if x in ["<","<=",">",">=","==","!="]]
	keys = [filter_string[i] for i in key_indis]
	return keys, key_indis

''' test
get_keys(["A", ">", "100", "and", "B", "==", "4"]) == (['A', 'B'], [0, 4])
'''


# extract information from info-field
def get_info_dict(info_field):
	'''
	in: INFO field from Vcf (string), e.g. "AC=27;AN=558;GF1=blabla"
	out: dictionary(key:value)
	'''
	info_field = info_field.split(";")
	info_dict = {}
	for info in info_field:
		i = info.split("=")
		# convert to number if possible
		try:
			info_dict[i[0]] = float(i[1])
		except ValueError:
			info_dict[i[0]] = i[1]
	return info_dict

''' test
get_info_dict("AC=27;AN=558;AF=0.048;GF0=0;GF1=0")
get_info_dict("AC=27;AN=558;AF=0.048;GF0=0;GF1=blabla")

'''


# use a string of expression to filter based on QUAL, FILTER and INFO column (needs whitespaces!)
def filter_with_string(filter_string, vcf_line):
	''' 
	in: filter_string (e.g. "QUAL > 0 and ( AN < 500 or AF == 0.048)") with mandatory whitespaces!
	    vcf_line as list of strings
	out: True/False
	'''
	filter_string = filter_string.split()
	keys, key_indis = get_keys(filter_string)
	# parse QUAL and FILTER column
	filter_dict = {"QUAL":float(vcf_line[5]), "FILTER":vcf_line[6]}
	# parse INFO-field only if necessary
	if not all([k in ["QUAL","FILTER"] for k in keys]):
		#print("Include filter on info field")
		filter_dict.update(get_info_dict(vcf_line[7]))
	# replace keywords with their filter_dict entry
	for i in key_indis:
		filter_string[i] = "filter_dict['" + filter_string[i] + "']"
	info_state = " ".join(filter_string)
	#print(filter_dict)
	#print(info_state)
	return eval(info_state)

''' test
vcf_line1 = ['21', '14896202', '.', 'C', 'T', '0', '.', 'AC=27;AN=558;AF=0.048;GF0=0;GF1=0', 'GT:GF']
vcf_line2 = ['21', '14896202', '.', 'C', 'T', '0', 'LowQual', 'AC=27;AN=558;AF=0.048;GF0=0;GF1=0', 'GT:GF']
vcf_line3 = ['21', '14896202', '.', 'C', 'T', '100', 'PASS', 'AC=27;AN=558;AF=a_string']

filter_with_string("QUAL > 0 and AC == 27", vcf_line1) == False
filter_with_string("( QUAL > 0 or FILTER != 'LowQual') and AC == 27", vcf_line1) == True
filter_with_string("( QUAL > 0 and FILTER == 'PASS')", vcf_line1) == False
filter_with_string("( QUAL > 0 or FILTER != 'LowQual') or ( AC == 27 and AF > 0.05)", vcf_line2) == False
filter_with_string("( QUAL > 0 or FILTER != 'LowQual') or ( AC == 27 and AF > 0.04)", vcf_line2) == True
filter_with_string("( AC <= 28 and AF >= 0.048)", vcf_line2) == True
filter_with_string("AF == 'bla' or ( QUAL > 0 and FILTER == 'PASS')", vcf_line3) == True
filter_with_string("AF == 'bla' and ( QUAL > 0 and FILTER == 'PASS')", vcf_line3) == False
filter_with_string("AF == 'a_string' and ( QUAL > 0 and FILTER == 'PASS')", vcf_line3) == True

'''



### filter based only on INFO-field
def filter_on_info(filter_string, info_dict):
	''' 
	in: filter_string (e.g. "AN < 500 or AF == 0.048") with mandatory whitespaces!
	    info_dict: dict[info_key:value]
	out: True/False
	'''
	filter_string = filter_string.split()
	# find keywords (positions before relational operator)
	key_indis = [i-1 for i, x in enumerate(filter_string) if x in ["<","<=",">",">+","==","!="]]
	for i in key_indis:
		filter_string[i] = "info_dict['" + filter_string[i] + "']"
	info_state = " ".join(filter_string)
	return eval(info_state)

''' test
info_dict = {'AF': 0.048, 'AC': 27.0, 'GF0': 0.0, 'AN': 558.0, 'GF1': "a_string"}
filter_on_info("AN < 500 and AF == 0.048", info_dict) == False
filter_on_info("AN < 500 or AF == 0.048", info_dict) == True
filter_on_info("AF > 0.04", info_dict) == True
filter_on_info("GF0 == 0 and AN > 550", info_dict) == True
filter_on_info("GF0 == 0.0 and AN > 550", info_dict) == True
filter_on_info("GF0 == 0.0 and AN > 550.0", info_dict) == True
filter_on_info("GF0 == 0 and AN > 560", info_dict) == False
filter_on_info("GF0 == 0 and GF1 == 'a_string'", info_dict) == True
filter_on_info("GF0 == 0 and GF1 != 'a_string'", info_dict) == False

'''





