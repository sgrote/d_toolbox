
'''
filter_line:
take a vcf-line and ONE string with a filter expression per line
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
def filter_line(filter_string, vcf_line):
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

filter_line("QUAL > 0 and AC == 27", vcf_line1) == False
filter_line("( QUAL > 0 or FILTER != 'LowQual') and AC == 27", vcf_line1) == True
filter_line("( QUAL > 0 and FILTER == 'PASS')", vcf_line1) == False
filter_line("( QUAL > 0 or FILTER != 'LowQual') or ( AC == 27 and AF > 0.05)", vcf_line2) == False
filter_line("( QUAL > 0 or FILTER != 'LowQual') or ( AC == 27 and AF > 0.04)", vcf_line2) == True
filter_line("( AC <= 28 and AF >= 0.048)", vcf_line2) == True
filter_line("AF == 'bla' or ( QUAL > 0 and FILTER == 'PASS')", vcf_line3) == True
filter_line("AF == 'bla' and ( QUAL > 0 and FILTER == 'PASS')", vcf_line3) == False
filter_line("AF == 'a_string' and ( QUAL > 0 and FILTER == 'PASS')", vcf_line3) == True

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


###############################

'''
filter_genotypes:
take a vcf-line and ONE string with a filter expression for SINGLE GENOTYPES per line
e.g. "GT=='0|0' or GF>=-1"
may contain all keywords described in FORMAT
return a modified vcf-line with non-passing genotypes replaced by "./." or ".|."

CAUTION: like in filter_line this needs whitespaces around keywords
(for easier and less error-prone spltting of expression and replacing the keywords)
(e.g. str.replace() can cause problems when a keyword is substring of another keyword)
'''


# replace keywords in filter_string with postion in format
# (for INFO-field above values are stored in a dictionary, here just index of key is used)
def get_condition_string(filter_string, gt_name):
	'''
	in: filter_string (with WHITESPACES),
	    gt_name: variable name for a gt_dict to evaluate expression on
	out: conditional expression to be used on a single genotype gt_name
	'''
	filter_string = filter_string.split()
	keys, key_indis = get_keys(filter_string)
	# replace keywords with their filter_dict entry
	for i in key_indis:
		filter_string[i] = gt_name + "['" + str(filter_string[i]) + "']"
	condition_string = " ".join(filter_string)
	return condition_string

''' test
get_condition_string("GF >= -1 or GT == '0|0'", "gt")
get_condition_string("( GF >= -1 or GT == '0|0') and SP == 5", "gt")
'''

# extract information from Genotypes given the format key values
def get_gt_dict(gt_parts, key_pos):
	'''
	in: list with entries from GT field, e.g. "[0/1", "-1"] from gt.split(":"); 
	    key_pos: list with entires from FORMAT-field
	out: dictionary(key:value)
	'''
	gt_dict = {}
	for i in range(len(gt_parts)):
		try:
			gt_dict[key_pos[i]] = float(gt_parts[i])
		except ValueError:
			gt_dict[key_pos[i]] = gt_parts[i]
	return gt_dict

''' test
get_gt_dict(["0/1","-1"], ["GT","GF"])

'''


# replace a vcf-line with filtered individual genotypes
def filter_genotypes(filter_string, vcf_line, gt_only=False):
	'''
	in: filter_string, e.g. GF >= 0 or GT == './0'
	    vcf_line as list of strings
	    if gt_only: leave only genotype, also adjust FORMAT to GT only
	out: vcf_line with individual genotypes not passing filter replaced with ./.
	'''
	condition_string = get_condition_string(filter_string, "gt_dict")
	out_line = vcf_line[:]
	key_pos = out_line[8].split(":")
	gt_pos = key_pos.index("GT")
	if gt_only:
		out_line[8] = "GT"
	# loop over fields with individual genotypes and evaluate expression
	for i in range(9, len(out_line)):
		gt_parts = out_line[i].split(":")
		gt_dict = get_gt_dict(gt_parts, key_pos)
		gt_pass = eval(condition_string)
		if not gt_pass:
			gt_parts[gt_pos] = "./."
			out_line[i] = ":".join(gt_parts)
		if gt_only:
			out_line[i] = gt_parts[gt_pos]
	return out_line

''' test
vcf_line = ['21','148','.','C','T','0','.','.','GT:GF','0/1:-1','0/1:3','./0:-1','.:-1','1/1:5','0/0:3']
filter_genotypes("GF >= 0 or GT == './0'", vcf_line) == ['21', '148', '.', 'C', 'T', '0', '.', '.', 'GT:GF', './.:-1', '0/1:3', './0:-1', './.:-1', '1/1:5', '0/0:3']
filter_genotypes("GF > 0", vcf_line, True) == ['21', '148', '.', 'C', 'T', '0', '.', '.', 'GT', './.', '0/1', './.', './.', '1/1', '0/0']

vcf_line2 = ['21','148','.','C','T','0','.','.','ST:GT:GF','P:0/1:-1','Q:0/1:3','P:./0:-1']
filter_genotypes("ST == 'P'", vcf_line2, True) == ['21', '148', '.', 'C', 'T', '0', '.', '.', 'GT', '0/1', './.', './0']
filter_genotypes("ST == 'Q' or GT == './0'", vcf_line2, True) == ['21', '148', '.', 'C', 'T', '0', '.', '.', 'GT', './.', '0/1', './0']
filter_genotypes("GF > 2 or GT == './0'", vcf_line2) == ['21', '148', '.', 'C', 'T', '0', '.', '.', 'ST:GT:GF', 'P:./.:-1', 'Q:0/1:3', 'P:./0:-1']
'''


# replace a vcf-line with genotypes only (without filtering)
def genotypes_only(vcf_line):
	'''
	in: vcf_line as list of strings
	out: vcf_line with genotypes only
	'''
	out_line = vcf_line[:]
	key_pos = out_line[8].split(":")
	gt_pos = key_pos.index("GT")
	out_line[8] = "GT"
	for i in range(9, len(out_line)):
		gt_parts = out_line[i].split(":")
		out_line[i] = gt_parts[gt_pos]
	return out_line

''' test
vcf_line = ['21','148','.','C','T','0','.','.','GT:GF','0/1:-1','0/1:3','./0:-1']
genotypes_only(vcf_line) == ['21', '148', '.', 'C', 'T', '0', '.', '.', 'GT', '0/1', '0/1', './0']
vcf_line2 = ['21','148','.','C','T','0','.','.','XX:GT:GF','x:0/1:-1','2/2:0/1:3']
genotypes_only(vcf_line2) == ['21', '148', '.', 'C', 'T', '0', '.', '.', 'GT', '0/1', '0/1']
vcf_line3 = ['21','148','.','C','T','0','.','.','GT','0/1','1|1']
genotypes_only(vcf_line3) == ['21', '148', '.', 'C', 'T', '0', '.', '.', 'GT', '0/1', '1|1']
'''


def check_filter(vcf_line, filter_string=None):
	''' True if line passes filter or no filter specified '''
	if filter_string and not filter_line(filter_string, vcf_line):
		return False
	else:
		return True


def filter_gt(vcf_line, filter_ind=None, gt_only=False):
	''' filter individual genotypes if specified, else return the input '''
	if filter_ind:
		vcf_line = filter_genotypes(filter_ind, vcf_line, gt_only) # this does gt_only on the fly
	elif gt_only:
		vcf_line = genotypes_only(vcf_line)
	return vcf_line[:]

''' test
vcf_line = ['21','148','.','C','T','0','.','.','GT:GF','0/1:-1','0/1:3','./0:-1']
filter_gt(vcf_line) == vcf_line
filter_gt(vcf_line, gt_only=True) == vcf_line[:8] + ['GT','0/1','0/1','./0']
filter_gt(vcf_line, filter_ind="GF >= 0", gt_only=True) == vcf_line[:8] + ['GT','./.','0/1','./.']
'''
