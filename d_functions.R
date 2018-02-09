
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
	limit = 8 # TODO: das wird hier noch nicht gebraucht da zu niedrig ==> aber code erstmal behalten
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



# plot barplot for D(pop1, pop2, pop3, pop4) with 3 pops fixed combi
plot_d_bars = function(input, superpops, ymin=NULL, ymax=NULL, mcex=0.9, legcex=0.8, auto_switch=TRUE){
	
	# switch pop1-pop2 if median D is negative
	if (auto_switch && (median(input$d) < 0)){
		input$d = input$d * -1
		input$z = input$z * -1
		input[,c(1,2)] = input[,c(2,1)]
	}
	
	# check which pop-position is variable 
	n_pops = apply(input[,1:4], 2, function(x) length(unique(x))) 
	varcol = which(n_pops != 1)
	if (length(varcol) != 1){
		stop("More than one of [pop1, pop2, pop3, pop4] is not unique.")
	}
	
	# title
	if (varcol == 1){
		titel = paste0("D(X, ",unique(input$pop2),", ",unique(input$pop3),", ",unique(input$pop4),")")
	} else if (varcol == 2){
		titel = paste0("D(",unique(input$pop1),", X, ",unique(input$pop3),", ",unique(input$pop4),")")
	} else if (varcol == 3){
		titel = paste0("D(",unique(input$pop1),", ",unique(input$pop2),", X, ",unique(input$pop4),")")
	} else if (varcol == 4){
		titel = paste0("D(",unique(input$pop1),", ",unique(input$pop2),", ",unique(input$pop3),", X)")
	}
	
	# subtitel
	subtitel = paste0(min(input$n_sites)," - ", max(input$n_sites), " informative sites")

	# convert to %
	input$se = 100 * (input$d / input$z)
	input$d = 100 * input$d
	
	# replace se=NA (caused by D=0)
	input$se[is.na(input$se)] = 0

	# merge with superpop-info and define spaces between bars
	input = cbind(input, superpops[match(input[,varcol], superpops[,1]), 2:4])
	input = input[order(input$order, input[,varcol]),]
	superbreak = c("dummy", input$super[-nrow(input)]) != input$super # for spaces
	spaces = rep(0.2, nrow(input))
	spaces[superbreak] = 0.8

	# dynamic y-axis
	if (is.null(ymin)) 	{
		ymax = max((max(input$d + input$se)), 0)
		ymin = min(min(input$d - input$se), 0)
		span = ymax - ymin
		ymax = ymax + abs(0.2*span)
		ymin = ymin - abs(0.05*span)
	}
	
	# cex xlab names (0.5-1)
	cexnames = 1 - (1/300) * nrow(input)
	cexnames = max(0.5, cexnames)

	# x-axis (this approach would create an additional plot in pdf since type="n" is not available for barplot) TODO: maybe change to plot(..., type="h")
	# and use this xlim in plot together with xaxs="i"
#	# do a fake plot to get the x-coords of the bars
#	bars = barplot(input$d, width=0.8, space=spaces, xaxs="i")
#	side_space = (bars[nrow(bars),1] - bars[1,1]) * 0.05
#	xlim=c(bars[1,1] - side_space, bars[nrow(bars),1] + side_space)
		
	input$asterix = ""
	input[abs(input$z) > 2, "asterix"] = "*"
	input[abs(input$z) > 3, "asterix"] = "**"
	
	# empty plot (to put horizontal lines before bars), panel.first didn't work with barplot
	barplot(input$d, col="white", border="white", ylim=c(ymin,ymax), width=0.8, space=spaces, xaxt="n", yaxt="n", xlab="")
	
	# horizontal lines
#	grid (NA, NULL, col="lightblue", lty=1) # NA: no lines on x-axis, NULL: autom lines at tick-marks on y
	yaxp = par("yaxp")
	hlines = seq(yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3])
	abline(h=hlines[-length(hlines)], col="lightblue", lty=1) # omit top line
	
	# barplot
	bars = barplot(input$d, names.arg=input[,varcol], col=colors()[input$color], main=titel, ylim=c(ymin,ymax), las=2, width=0.8, space=spaces, ylab="D[%]", cex.lab=mcex, cex.axis=0.8, cex.names=cexnames, cex.main=1.2,  add=T)
#	# xlab
#	mtext("X", 1, 5.5, cex=mcex)

	# error bars
	suppressWarnings(arrows(c(bars,bars), c(input$d+input$se,input$d-input$se), c(bars,bars), c(input$d, input$d), angle=90, code=1, length=0.015))
	# significance asterix
	input$bars = bars
	posi = input[input$d >= 0,]
	nega = input[input$d < 0,]
	if(nrow(posi) > 0){
		text(x=posi$bars, y=posi$d+posi$se, labels=posi$asterix, pos=3, offset=0, cex=0.5)
	}
	if(nrow(nega) > 0){
		text(x=nega$bars, y=nega$d-nega$se, labels=nega$asterix, pos=1, cex=0.5)		
	}
	
	# subtitle
	mtext(subtitel, line=-1, adj=0.5, cex=mcex) 
#	# label on y-axis ends # TODO: not applicable anymore now that also pop1 or pop2 could be X-axis
#	mtext(substring(pop1,1,3), side=3, adj=-0.03, cex=mcex)
#	mtext(substring(pop2,1,3), side=1, adj=-0.03, cex=mcex)

	# legend Z-scores
	legend("topleft",cex=legcex,ncol=2,legend=c("*", "**", "|Z| > 2","|Z| > 3"),bty="n",title="weighted block Jackknife")
	# legend populations
	legend("topright",cex=legcex,ncol=3,legend=unique(input$super),bty="n",fill=colors()[unique(input$color)],title="super populations")

}


