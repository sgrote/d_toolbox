#!/usr/bin/env Rscript

# compute ABBA and BABA counts for D(pop1, pop2, pop3, pop4) per ANCESTRAL-DERIVED allele pair
# based on RANDOMLY drawn alleles 
# (TODO: maybe draw a random allele only for the outgroup,
#    or make it accurate with A/B switching for non-binary outgroup, but probably it doesnt matter)
# needs 'out_d' and a full 'sites' file (run d_stats.py with --sites 'full' and d_genomewide.R before)
# population-matches are retrieved from out_d, can also take another file which might be a subset pop-matches
# output: abba,baba,per pop-match and A/B allele (A,C,T,G)

library(optparse)
library(readr)
library(gtools)

option_list = list(
	make_option(c("-d", "--infile"), type="character", default="out_d",
		help="d-stats input \n\t\tdefault = %default"),
	make_option(c("-o", "--outfile"), type="character", default="abba_baba_per_base",
		help="output file name \n\t\tdefault = %default"),
	make_option(c("-t", "--transver"), action="store_true", default=FALSE,
		help="input has transversions only")
)

opt_parser = OptionParser(option_list=option_list, description="\ncompute ABBA and BABA counts for D(pop1, pop2, pop3, pop4) per ancestral/derived (A/B) allele pair
needs 'out_d' and a full 'sites' file (run d_stats.py with --sites 'full' and d_genomewide.R before)
population-matches are retrieved from out_d, can also take another file which might be a subset of pop-matches
output: abba, baba per pop-match and A/B allele (A,C,T,G)")
opt = parse_args(opt_parser)
print(opt)


####### helper

# pop1 pop2 pop3 pop4 A B abba baba  per pop1-pop4 combi and A_B allele pattern
count_abba_baba = function(freqs, pops, base_combis){
	
	# subset to pop-combi
	freqqi = cbind(freqs[,2:4], freqs[,pops])
	# remove NA
	freqqi = freqqi[!is.na(rowSums(freqqi[,pops])),]
	# fix pop4 as A
	pop4_ALT = freqqi[, pops[4]] == 1
	freqqi[pop4_ALT, pops] = 1 - freqqi[pop4_ALT, pops]
	freqqi[pop4_ALT, c("ref","alt")] = freqqi[pop4_ALT, c("alt","ref")]
	# count ABBA/BABA
	abba = freqqi[,pops[1]] == 0 & freqqi[,pops[2]] == 1 & freqqi[,pops[3]] == 1
	baba = freqqi[,pops[1]] == 1 & freqqi[,pops[2]] == 0 & freqqi[,pops[3]] == 1
	# split by A/B
	a_b = paste(freqqi$ref, freqqi$alt, sep="_")
	abba_sum = tapply(abba, a_b, sum)
	baba_sum = tapply(baba, a_b, sum)
	# put in the right order, replace NA by 0,
	abba_sum = abba_sum[match(base_combis, names(abba_sum))]
	baba_sum = baba_sum[match(base_combis, names(baba_sum))]
	county = cbind(abba_sum, baba_sum)
	county[is.na(county)] = 0
	row.names(county) = base_combis
	
	return(county)
}

####### main

# get all base-combinations
bases = c("A","C","T","G")
bc = data.frame(A_allele=bases[rep(1:4, each=4)], B_allele=bases[rep(1:4, length.out=4)])
bc = bc[bc[,1] != bc[,2],]
bc$A_B = paste(bc[,1], bc[,2], sep="_")
if (opt$transver){
	bc = bc[!(bc$A_B %in% c("G_A","A_G","C_T","T_C")),]
}

# read out_d
if (! opt$infile %in% dir()){
	stop(opt$infile, " not found. Needs 'out_d'-like file created by d_stats.py")
}    
out_d = read.table(opt$infile, header=T, as.is=T)

# get population matches pop1-pop2-pop4
fixed_pops = unique(out_d[,1:4])

# initialize output
out = data.frame(fixed_pops[rep(1:nrow(fixed_pops), each=nrow(bc)),])
out = cbind(out, bc[,1:2]) # get recycled
out$abba = 0
out$baba = 0
row.names(out) = 1:nrow(out)

# pop-match strings
paste_pops = paste(out[,1], out[,2], out[,3], out[,4], sep="_")

# get chrom subdirs
subdirs = dir()[file.info(dir())$isdir]
chroms = subdirs[substr(subdirs,1,3)=="chr"]
chroms = mixedsort(chroms)

## loop over chroms
for (i in 1:length(chroms)){
	setwd(chroms[i])
	print(getwd())
	# check that sites-file is there and has right format, skip chrom with message if not
	if (! "sites" %in% dir()){
		message("skipping ", chroms[i], " (no sites file).")
		setwd("../")
		next
	}
	sites = suppressMessages(as.data.frame(read_tsv("sites")))
	if (ncol(sites) < 6){
		message("skipping ", chroms[i], ". sites file lacks columns (run d_stats.py with --sites 'full' option).")
		setwd("../")
		next
	}
	# convert to numeric (None -> NA)
	sites[,5:ncol(sites)] = suppressWarnings(sapply(sites[,5:ncol(sites)], as.numeric))
	# draw random allele
	sites[,5:ncol(sites)] = apply(sites[,5:ncol(sites)], 2, function(x){
		suppressWarnings(rbinom(length(x), 1, x))})

	## loop over pops
	for (p in unique(paste_pops)){
		pops = unlist(strsplit(p, "_"))
		# pop1 pop2 pop3 pop4 A B abba baba  per pop1-pop4 combi
		counts = count_abba_baba(sites, pops, bc[,3])
		out[paste_pops == p, c("abba","baba")] = out[paste_pops == p, c("abba","baba")] + counts
	}

	setwd("../")
}

write.table(out, opt$outfile, row.names=F, quote=F, sep="\t")





























