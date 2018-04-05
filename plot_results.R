
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


# plot d in heatmap given dataframe: [pop1, pop2, pop3, pop4, d, z] for one pop3
plot_d = function(dtab, sub="", dendro=FALSE){
	if(length(unique(dtab[,3])) != 1 | length(unique(dtab[,4])) != 1){
		stop("Multiple populations in column 3 or 4")
	}
	# convert to D[%]
	dtab[,5] = 100*dtab[,5]
	# initialize d-matrix
	pops = unique(c(as.matrix(dtab[,1:2])))
	sqd = matrix(nrow=length(pops),ncol=length(pops))
	colnames(sqd) = pops
	rownames(sqd) = pops
	# initialize z-matrix
	sqz = sqd
	# populate with d- and z-values
	for (i in 1:nrow(dtab)){
		sqd[dtab[i,1], dtab[i,2]] = dtab[i,5] 
		sqd[dtab[i,2], dtab[i,1]] = -dtab[i,5] 
		sqz[dtab[i,1], dtab[i,2]] = dtab[i,6] 
		sqz[dtab[i,2], dtab[i,1]] = -dtab[i,6] 
	}		
	# replace NA with 0 (diagonal)
	sqd[is.na(sqd)] = 0
	sqz[is.na(sqz)]= 0
	# NEW: significance asterix-matrix
	sqa = matrix(ncol=ncol(sqz), nrow=nrow(sqz), data="", dimnames=list(rownames(sqz),colnames(sqz)))
#	sqa[abs(sqz) > 2] = "*"
#	sqa[abs(sqz) > 3] = "**"
	sqa[sqz > 2] = "*"
	sqa[sqz	 > 3] = "**"
	# add custom breaks to allow fixed limits 
	limit = 8 # TODO: das wird hier noch nicht gebraucht da zu niedrig ==> aber code erstmal behalten
	if (max(sqd) > limit) limit = max(sqd)
	breakys = seq(-limit, limit, len=101)
	# titel und so
	titel = paste("D(Row, Col, ", dtab[1,3],", ",dtab[1,4],")", sep="")
	cexRow = (0.2 + 1/log10(nrow(sqd))) * 0.7
	dtype = ifelse(dendro, "row", "none") # row and column dendrogramm are symmetric
	gplots::heatmap.2(sqd, main=titel,
		scale="none", density.info="none", trace="none",
		cexRow=cexRow, cexCol=cexRow, keysize=1.2, margins=c(6,6),
		breaks=breakys, col=heat.colors(100), symkey=F,
		dendrogram=dtype, Rowv=dendro, Colv=dendro,		 
		cellnote=sqa, notecol="black", notecex=1.2) # asterix
	# legend z-score mark
	legend("top",cex=0.7,ncol=2,legend=c("*","**","|Z| > 2","|Z| > 3"), bty="n",pt.bg=1,title="weighted block Jackknife")
#	legend(x=0.85,y=1.1,cex=0.7,ncol=2,legend=c("*","**","|Z| > 2","|Z| > 3"), bty="n",pt.bg=1,title="weighted block Jackknife")
	# subtitle
	legend("topright", sub, bty="n")   	
}



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

