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


def get_pops_from_files(pop1, pop2, pop3, pop4):
    '''
    in: files with one column containing population names
    out: list [[pop1][pop2][pop3][pop4]] for all matches,
     avoiding self-matches and duplicates (including D(w,x,y,z), D(x,w,y,z))
    '''
    # TODO: maybe use a pandas df here?
    pop_files = [pop1, pop2, pop3, pop4]
    pops = [[] for k in range(4)]
    pw_pops = [[] for k in range(4)]
    # read populations
    for i in range(4):
        with open(pop_files[i], "r") as p:
            for pop in p:
                pops[i].append(pop.rstrip())
    print(pops)
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
    if len(pw_pops[1]) == 0:
        sys.exit("\nError: could not create pairwise population matches - check the pop input files.")
    return pw_pops


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
    

# NEW use pandas here instead of get_pop_dict(), allow multiple pops per sample
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
    for i in range(len(vcf_header)):
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
            # TODO maybe code gender as 0/1, m/f
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

# TODO: make this faster
def get_p(vcf_line, pop, col_gender, X=False):
    '''
    alternative allele-freqs for one population for one biallelic site
    in: vcf = vcf-line as list
        pop = [colnumbers] for the population
        col_gender = dictionary {colnumber:sex}
        X = T/F if its X-chrom
    out: ALT-allele-freq
    '''
    ac = 0
    sum_alleles = 0
    # count alternative alleles and total number of alleles for one population  
    for colnumber in pop:
        # X-ch for males (handle 1|., 0, ., ./., 1|1, and so on) (take X as input to not check every line)
        if X and col_gender[colnumber]=="male":
            # alternative alleles
            if "1" in vcf_line[colnumber]:
                ac += 1
            # all alleles
            #if not(vcf_line[colnumber] == "." or vcf_line[colnumber] == ".|."):
            if not all([char in {"/",".","|"} for char in vcf_line[colnumber]]):
                sum_alleles += 1
        else:
            # alternative alleles
            if vcf_line[colnumber] in {"1|1", "1/1"}:
                ac += 2
            elif "1" in vcf_line[colnumber]:
                ac += 1
            # all alleles
            if "." not in vcf_line[colnumber]:
                sum_alleles += 2
            elif vcf_line[colnumber] not in {".|.", "./."}:
                sum_alleles += 1
    # compute frequency 
    if sum_alleles == 0:
        #print("This Pop has no genotypes:", pop)
        return None   # unknown sites may be fully genotyped invariable sites, here there is simply no GT 
    else: 
        daf = ac*1.0 / sum_alleles
    return daf

''' test
vcf_line1 = ['21','148','.','C','T','0','.','.','GT','0/1','1/1','./0','./.','1|.','1|1','.|.']
vcf_line2 = ['X','148','.','C','T','0','.','.','GT','1/1','1','./0','./.','1|.','1/0','.']
col_gender = {9:"male", 10:"male", 11:"female", 12:"female", 13:"male", 14:"female", 15:"male"}

get_p(vcf_line1, [9,10], col_gender) == 0.75
get_p(vcf_line1, [9,10,11], col_gender) == 3.0/5
get_p(vcf_line1, [9,10,11,14,15], col_gender) == 5.0/7
get_p(vcf_line1, [12,13,14,15], col_gender) == 1
get_p(vcf_line2, [9,10], col_gender, True) == 1
get_p(vcf_line2, [9,11], col_gender, True) == 0.5    # count male only once for X even if 1/1
get_p(vcf_line2, [9,14], col_gender, True) == 2.0/3  # but still count female twice
get_p(vcf_line2, [9,10,11], col_gender, True) == 2.0/3
get_p(vcf_line2, [9,10,11,12], col_gender, True) == 2.0/3

'''

    

# one list for containig block-abba-baba-sums for altai and vindija together in one list
# 0. blocks, 1. abba, baba, 2. pw_pops  [ [[abba][baba]] [[abba][baba]] ...]
# CUATION this assumes pure genotype-input
def abba_block_sums(vcf, pop_colnums, pw_pops, col_gender, block_size, centro_range=None, transver=False, sitesfile=None):
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
    # initialize abba and baba block-sums [[abba],[baba]] 
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
        ### TODO: remove
        #if int(line[1]) > 17419425:
            #continue
        #print(line)
        # get first pos for block-border and check if it's chrX
        if first:
            block_end = int(line[1]) + block_size
            X = True if line[0]=="X" else False
            first = False
        # exclude centromere (usually not needed because of manifesto)
        if centro_range and (int(line[1]) > centro_range[0] and int(line[1]) < centro_range[1]):
            #print("Skipping centromere position", int(line[1]), centro_range)
            continue
        ## bases/genotypes
        ref = line[3]
        alt = line[4]
        ## reduce to biallelic 
        # TODO: add fancy biallelic version that evaluates biallelic for pops of interest or even comps
        if alt not in "ACTG":
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
            block_sum = [[0]*n_comp, [0]*n_comp]
            block_end = block_end + block_size
            blocksites = 0
            block_count += 1
            print("starting block ", block_count+1, "until position", block_end)

        ## TODO: maybe before computing p(ALT), check that not all (pop3 and pop4) or (pop1 and pop2) genotypes are the same
        ## TODO: maybe also check if in the genotypes of e.g. pops1 or pops4 any "1" or "0" is present
        ## get alternative allele freqs for input populations
        p = {}
        for pop in pops:    
            p[pop] = get_p(line, pop_colnums[pop], col_gender, X)
        #print(p)
        ## check that not all pop-freqs are 0 or 1 and None (pop with variant might not be part of D-stats)
        pset = set(p.values())
        if pset.issubset({0,None}) or pset.issubset({1,None}):
            #print("Uninformative site:", pset) 
            continue

        ## compute p(ABBA) and p(BABA)  
        abba = [] # TODO: initialize with n_comp and 0 and then acces with index instead of append and just skip if None in ....
        baba = []
        for i in range(n_comp):
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
                # add to sites-count-per pw_pop if informative
                if p_abba > 0 or p_baba > 0:
                    #print(i)
                    #print(sites_comp)
                    sites_comp[i] += 1
        # check that not all is 0
        if sum(abba)==0 and sum(baba)==0:
            #print("uninformative site", line)
            continue
        # add abba and baba to block_sum
        for i in range(n_comp):
            block_sum[0][i] += abba[i]
            block_sum[1][i] += baba[i]  
        ## add to sites counter and TODO sites file
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
