
# go through all autosome subdirectories and combine abba and baba 
# genomewide only uses autosomes!

# jackknife:
source("/mnt/expressions/steffi/D/d_toolbox/d_functions.R")
require(gtools)
require(plyr)

# get chrom subdirs
subdirs = dir()[file.info(dir())$isdir]
chroms = subdirs[substr(subdirs,1,3)=="chr"]
chroms = mixedsort(chroms)

# initialize genomewide
genome_blocks = data.frame()
# go through directories
for (i in 1:length(chroms)){
	setwd(chroms[i])
	print(getwd())
	## TODO: only use sites file if one is present
	sites = as.integer(strsplit(system("wc -l sites", intern=T), " ")[[1]][1])	
	out_blocks = read.table("out_blocks", header=T)
	pw_sites = read.table("sites_comp", as.is=T, header=T)
	# jackknife - by-function also returned NULL for unused factor combinations 
	d = ddply(out_blocks, .(pop1,pop2,pop3,pop4), d_jackknife)
	# NEW: add sites per pw
	out = merge(d, pw_sites)
	write.table(out, "out_d", row.names=F, quote=F, sep="\t")	
	# collect for genomewide
	if (chroms[i] != "chrX"){
		if (i==1){
			genome_pw_sites = pw_sites
			n_sites_genome = data.frame("chr"=i, "sites"=sites)
			last_block = 0
		} else {
			genome_pw_sites$n_sites = genome_pw_sites$n_sites + pw_sites$n_sites
			last_block = genome_blocks[nrow(genome_blocks),"block"]
		}
		out_blocks$block = out_blocks$block + last_block	
		genome_blocks = rbind(genome_blocks, out_blocks)
		if (i!=1) {
			n_sites_genome = rbind(n_sites_genome, c(i, sites))
		}
	} else {
		message("skipping chrX for genomewide results")
	}	
	setwd("../")
}
# write genome-wide sites-files
write.table(n_sites_genome, "n_sites_genome", quote=FALSE, row.names=FALSE)
# genomewide D
genome_d = ddply(genome_blocks, .(pop1,pop2,pop3,pop4), d_jackknife) 
out = merge(genome_d, genome_pw_sites)
write.table(out, "out_d", row.names=F, quote=F, sep="\t")
