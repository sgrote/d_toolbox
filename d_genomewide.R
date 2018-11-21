#!/usr/bin/env Rscript


# go through all autosome subdirectories and combine abba and baba 
# genomewide only uses autosomes!

require(gtools)
require(plyr)

# jackknife:
argv = commandArgs(trailingOnly = FALSE)
base_dir = dirname(substring(argv[grep("--file=", argv)], 8))
source(file.path(base_dir, "d_functions.R"))

# get chrom subdirs
subdirs = dir()[file.info(dir())$isdir]
chroms = subdirs[substr(subdirs,1,3)=="chr"]
chroms = mixedsort(chroms)

# initialize genomewide
genome_blocks = data.frame()
# go through directories
first = TRUE
for (i in 1:length(chroms)){
	setwd(chroms[i])
	print(getwd())
	if (! "out_blocks" %in% dir()){
		message("skipping ", chroms[i], " (no out_blocks file).")
		setwd("../")
		next
	}
	sitesfile = FALSE
	if ("sites" %in% dir()){
		sitesfile = TRUE
		sites = as.integer(strsplit(system("wc -l sites", intern=T), " ")[[1]][1])
	}    
	out_blocks = read.table("out_blocks", header=T)
	pw_sites = read.table("sites_comp", as.is=T, header=T)
	# jackknife - by-function also returned NULL for unused factor combinations 
	d = ddply(out_blocks, .(pop1,pop2,pop3,pop4), d_jackknife)
	# NEW: add sites per pw
	out = merge(d, pw_sites)
	write.table(out, "out_d", row.names=F, quote=F, sep="\t")	
	# collect for genomewide
	if (chroms[i] != "chrX"){
		if (first){
			if (sitesfile){
				n_sites_genome = data.frame("chr"=i, "sites"=sites)
			}
			genome_pw_sites = pw_sites
			last_block = 0
		} else {
			if (sitesfile){
				n_sites_genome = rbind(n_sites_genome, c(i, sites))
			}
			genome_pw_sites$n_sites = genome_pw_sites$n_sites + pw_sites$n_sites
			last_block = genome_blocks[nrow(genome_blocks),"block"]
		}
		out_blocks$block = out_blocks$block + last_block	
		genome_blocks = rbind(genome_blocks, out_blocks)
		first = FALSE
	} else {
		message("skipping chrX for genomewide results")
	}	
	setwd("../")
}
# write genome-wide sites-files
if (sitesfile){
	write.table(n_sites_genome, "n_sites_genome", quote=FALSE, row.names=FALSE)
}
# genomewide D
genome_d = ddply(genome_blocks, .(pop1,pop2,pop3,pop4), d_jackknife) 
out = merge(genome_d, genome_pw_sites)
write.table(out, "out_d", row.names=F, quote=F, sep="\t")

