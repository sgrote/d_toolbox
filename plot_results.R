
# make a heatmap for of D(X, Y , pop3, pop4) for every pop3-pop4 combinations

# Plot D(pop1, pop2, X, pop4) for every pop1-pop2-pop4 combi as barplot

### command line arguments:
# 1) infile (out_d)
# 2) outfile.pdf

args = commandArgs(trailingOnly=TRUE)
out_d_file = args[1]
out_pdf = args[2]


# TODO: instead of loop over pops check present combinations (maybe  not all pairwise are present)
# TODO: mover other plotting scripts from ../scripts to ../old_scripts

source("/mnt/expressions/steffi/D/d_toolbox/d_functions.R") # plot_d


######################################
out_d = read.table(out_d_file, as.is=T, header=T)

## plot genomewide results
pdf(out_pdf)
	par(cex.main=0.8, oma=c(1,0,0,1.5), cex.lab=1)
	combis = unique(out_d[,c("pop3", "pop4")])
	for (i in 1:nrow(combis)){
		one34 = out_d[out_d$pop3==combis[i,"pop3"] & out_d$pop4==combis[i,"pop4",],]
		# heatmap
		plot_d(one34, dendro=TRUE)
	}
dev.off()

