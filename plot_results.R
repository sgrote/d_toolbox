
# D_heatmaps.pdf
# make a scatterplot for of D(X, Y , pop3, pop4) for every pop3-pop4 combinations

source("/mnt/expressions/steffi/D/d_toolbox/d_functions.R") # plot_d


######################################
out_d = read.table("out_d", as.is=T, header=T)

## TODO: add a pop4 loop, since there might be rare cases with several outgroups in one analysis
## plot genomewide results
pdf("D_heatmaps.pdf")
	par(cex.main=0.8, oma=c(0.5,0,0,1), cex.lab=1)
	for (pop3 in unique(out_d$pop3)){
		one3 = out_d[out_d$pop3==pop3,]
		# scatterplots
		plot_d(one3, dendro=TRUE)
	}
dev.off()

