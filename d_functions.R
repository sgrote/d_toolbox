
library(gplots)

# compute jackknife given all blocks, one pop1-pop2-pop3-pop4 combination
d_jackknife = function(block_data){
	
	# remove uninformative blocks
	block_data = block_data[block_data[,"abba"] + block_data[,"baba"] != 0,]
	# account for 'all blocks empty'
	if(nrow(block_data)==0){
		# return D and Z-score
		out = c("d"=0, "z"=0, "abba"=0, "baba"=0)	
		return(out)
	}	
	# variables like in Thesis/Busing_1999
	g = nrow(block_data) 
	n = sum(block_data[,c("abba","baba")])
	mj = block_data[,"abba"] + block_data[,"baba"]
	h = n / mj

	# compute D for all blocks together (D_hat)
	full_abba = sum(block_data[,"abba"])
	full_baba = sum(block_data[,"baba"])
	D_hat_full = (full_baba - full_abba) / (full_baba + full_abba)

	# compute D leaving out one block at a time (D_hat_j_star) and weight
	loo_abba = full_abba - block_data[,"abba"]
	loo_baba = full_baba - block_data[,"baba"]
	D_hat_j_star = (loo_baba - loo_abba) / (loo_baba + loo_abba)
	# jackknife estimator
	D_hat_J = g * D_hat_full - sum((1-mj/n) * D_hat_j_star)
	# pseudo-observations D_wave
	D_wave = h * D_hat_full - (h-1) * D_hat_j_star
	
	# variance estimator
	vari = sum(1/(h-1) * (D_wave - D_hat_J)^2) / g
	# Z
	z = D_hat_full / sqrt(vari)
	# NEW: abba and baba counts
	# return D and Z-score
	out = c("d"=D_hat_full, "z"=z, "abba"=full_abba, "baba"=full_baba)
	return(out)
}



# plot sites in heatmap given dataframe: [pop1, pop2, pop3, pop4, sites] for one pop3
plot_sites = function(dtab, sub=""){
	if(length(unique(dtab[,3])) != 1 | length(unique(dtab[,4])) != 1){
		stop("Multiple populations in column 3 or 4")
	}
	# initialize d-matrix
	pops = unique(c(as.matrix(dtab[,1:2])))
	sqd = matrix(nrow=length(pops),ncol=length(pops))
	colnames(sqd) = pops
	rownames(sqd) = pops
	# populate with d- and z-values
	for (i in 1:nrow(dtab)){
		sqd[dtab[i,1], dtab[i,2]] = dtab[i,5] 
		sqd[dtab[i,2], dtab[i,1]] = dtab[i,5] 
	}
	# replace NA with 0 (diagonal)
	sqd[is.na(sqd)] = 0		
	# add custom breaks to allow fixed limits 
	limit = 8 # das wird hier noch nicht gebraucht da zu niedrig ==> aber code erstmal behalten
	if (max(sqd) > limit) limit = max(sqd)
	breakys = seq(0, limit, len=101)
	# titel und so
	titel = paste("sites: D(Row, Col, ", dtab[1,3],", ",dtab[1,4],")", sep="")
	cexRow = (0.2 + 1/log10(nrow(sqd))) * 0.7
	gplots::heatmap.2(sqd, main=titel,
		scale="none", density.info="none", trace="none",
		cexRow=cexRow, cexCol=cexRow, keysize=1.2, margins=c(6,6),
		breaks=breakys, col=heat.colors(100), symkey=F,
		Rowv=FALSE, Colv=FALSE, dendrogram="none")		 
	legend("topright", sub, bty="n")   	
}





