#!/usr/bin/env Rscript

# compute D(pop1, pop2, pop3, pop4) per derived freq-bin in pop3 (per pop1 or pop2 makes no sense)
# needs 'out_d' and a full 'sites' file (run d_stats.py with --sites 'full' and d_genomewide.R before)
# population-matches are retrieved from out_d, can also take another file which might be a subset pop-matches
# output: abba,baba,n_sites,d,z per pop-match and freq-bin (bin = mean of upper and lower bin-limit, bin=1 -> fixed)

library(optparse)
library(readr)
library(gtools)
library(plyr)


option_list = list(
	make_option(c("-d", "--infile"), type="character", default="out_d",
		help="d-stats input \n\t\tdefault = %default"),
	make_option(c("-o", "--outfile"), type="character", default="out_d_freqbins",
		help="output file name \n\t\tdefault = %default"),
	make_option(c("-b", "--bins"), type="integer", default=50,
		help="nr. of frequency bins. 'Fixed derived' is always the additional bin 1.0 \n\t\tdefault = %default")
)

opt_parser = OptionParser(option_list=option_list, description="\ncompute D(pop1, pop2, pop3, pop4) per derived freq-bin in pop3 (by pop1 or pop2 makes no sense)
needs 'out_d' and a full 'sites' file (run d_stats.py with --sites 'full' and d_genomewide.R before)
population-matches are retrieved from out_d, can also take another file which might be a subset of pop-matches
output: abba, baba, n_sites, d per pop-match and freq-bin (bin = mean of upper (exclusive) and lower (inclusive) bin-limit, bin=1 -> fixed)")
opt = parse_args(opt_parser)
print(opt)

# jackknife:
argv = commandArgs(trailingOnly = FALSE)
base_dir = dirname(substring(argv[grep("--file=", argv)], 8))
source(file.path(base_dir, "d_functions.R"))


## main

# create bins
bin_borders = seq(0, 1, length.out=opt$bins+1)
bin_width = bin_borders[2] - bin_borders[1]
# means of bin-borderes + fixed
bins = c(bin_borders[bin_borders != 1] + (bin_width/2), 1)

# read out_d
if (! opt$infile %in% dir()){
	stop(opt$infile, " not found. Needs 'out_d'-like file created by d_stats.py")
}    
out_d = read.table(opt$infile, header=T, as.is=T)

# get population matches pop1-pop2-pop4
out_d$fixed_pop = paste(out_d[,1], out_d[,2], out_d[,4], sep="_")
fixed_pops = unique(out_d$fixed_pop)

# initialize sums per freq-bin and pop1-pop2-pop3-pop4 combi (list)
out_d$paste_pops = paste(out_d[,1], out_d[,2], out_d[,3], out_d[,4], sep="_")
freq_list = vector("list", length=nrow(out_d))
names(freq_list) = out_d$paste_pops

# get chrom subdirs
subdirs = dir()[file.info(dir())$isdir]
chroms = subdirs[substr(subdirs,1,3)=="chr"]
chroms = mixedsort(chroms)

## loop over chrom-subdirectories
last_block = 0
for (i in 1:length(chroms)){
	setwd(chroms[i])
	print(getwd())
	# check that sites-file is there and has right format, skip chrom with message if not
	if (! "sites" %in% dir()){
		message("skipping ", chroms[i], " (no sites file).")
		setwd("../")
		next
	}
	sites = suppressMessages(as.data.frame(read_tsv("sites", col_types=cols(.default="c"))))
	if (ncol(sites) < 6){
		message("skipping ", chroms[i], ". sites file lacks columns (run d_stats.py with --sites 'full' option).")
		setwd("../")
		next
	}
	# convert blocks and freqs to numeric (None -> NA)
	sites[,c(1,5:ncol(sites))] = suppressWarnings(sapply(sites[,c(1,5:ncol(sites))], as.numeric))
	# increment block nrs, independent 1:x per chrom
	sites$block = sites$block + last_block
	last_block = sites[nrow(sites), "block"]
	## loop over pop1-pop2-pop4 population
	for (fp in fixed_pops){
		one124 = out_d[out_d$fixed_pop == fp,]
		pop1 = one124[1,1]
		pop2 = one124[1,2]
		pop4 = one124[1,4]
		# reduce outgroup to {0,1}
		freqs = sites[sites[,pop4] %in% c(0,1),]
		# remove NA in pop1/pop2
		freqs = freqs[! (is.na(freqs[,pop1]) | is.na(freqs[,pop2])),]
		# convert to derived allele freq (based on pop4)
		freqs[freqs[,pop4]==1, 5:ncol(freqs)] = 1-freqs[freqs[,pop4]==1, 5:ncol(freqs)]
		# compute p(AB_A) and p(BA_A) per sites
		p_ab = (1-freqs[,pop1]) * freqs[,pop2]
		p_ba = freqs[,pop1] * (1-freqs[,pop2])
		## for every pop3 get ABBA/BABA-sums per freq-bin
		for (j in 1:nrow(one124)){
			pop3 = one124[j,"pop3"]
			pp = one124[j,"paste_pops"]
			# remove pop3 fixed A or NA)
			b0 = freqs[,pop3] %in% c(0,NA)
			abba = p_ab[!b0] * freqs[!b0,pop3]
			baba = p_ba[!b0] * freqs[!b0,pop3]
			inform_site = (abba != 0 | baba != 0) * 1
			# get freq bins and blocks, bin = mean of borders (1 only for fixed)
			freqbins = floor(freqs[!b0,pop3]*opt$bins) / opt$bins
			freqbins[freqbins != 1] = freqbins[freqbins != 1] + (bin_width/2)
			blocks = freqs[!b0, "block"]
			# ABBA/BABA/n_sites per bin and block
			counts = aggregate(data.frame(abba, baba, inform_site), list(bin=freqbins, block=blocks), sum)
			# add to genomewide
			freq_list[[pp]] = append(freq_list[[pp]], list(counts))
		}
	}
	setwd("../")
}

# for each list element: combine chromosomes
freq_list = lapply(freq_list, function(x){do.call(rbind, x)})

# for each list element: run jackknife per bin
d = lapply(freq_list, function(x){ddply(x, .(bin), d_jackknife)})
n = lapply(freq_list, function(x){ddply(x, .(bin), function(y) {sum(y$inform_site)})})

# create output
out = cbind(do.call(rbind,d), do.call(rbind,n))
out = cbind(out_d[rep(1:nrow(out_d), each=nrow(d[[1]])),1:4], out)
out = out[,c(1:5,8,9,11,6,7)]
colnames(out)[8] = "n_sites"

## account for ABBA+BABA = 0
out$d[is.na(out$d)] = 0
out$z[is.na(out$z)] = 0

## save
write.table(out, opt$outfile, quote=F, row.names=F, sep="\t")



