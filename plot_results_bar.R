#!/usr/bin/env Rscript

# Plot D(pop1, pop2, X, pop4) for every pop1-pop2-pop4 combi as barplot
# group, order and color bars by info given in .csv file 

library("optparse")

#### helper

# plot barplot for D(pop1, pop2, pop3, pop4) with 3 pops fixed combi
plot_d_bars = function(input, superpops, ymin=NULL, ymax=NULL, mcex=0.9, legcex=0.8, lege=c(T,T), infsites=TRUE, namesfile=NULL){
	

	# check which pop-position is variable 
	n_pops = apply(input[,1:4], 2, function(x) length(unique(x))) 
	varcol = which(n_pops != 1)
	if (length(varcol) != 1){
		stop("More than one of [pop1, pop2, pop3, pop4] is not unique.")
	}

	# merge with superpop-info and define spaces between bars
	input = cbind(input, superpops[match(input[,varcol], superpops[,1]), 2:4])
	input = input[order(input$order, input[,varcol]),]
	superbreak = c("dummy", input$super[-nrow(input)]) != input$super # for spaces
	spaces = rep(0.2, nrow(input))
	spaces[superbreak] = 0.8

	# convert to official names
	if (! is.null(namesfile)){
		offinames = read.csv(namesfile, header=T, as.is=T)
		input[,1:4] = apply(input[,1:4], 2, get_offi, offinames)
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
	
	# convert to %
	input$se = 100 * (input$d / input$z)
	input$d = 100 * input$d
	
	# replace se=NA (caused by D=0)
	input$se[is.na(input$se)] = 0

	# dynamic y-axis
	if (is.null(ymin)) 	{
		ymax = max((max(input$d + input$se)), 0)
		ymin = min(min(input$d - input$se), 0)
		span = ymax - ymin
		ymax = ymax + abs(0.5*span)
		ymin = ymin - abs(0.2*span)
	}
	
	# cex xlab names (0.5-1)
	cexnames = 1 - (1/300) * nrow(input)
	cexnames = max(0.5, cexnames)

	input$asterix = ""
	input[abs(input$z) > 2, "asterix"] = "*"
	input[abs(input$z) > 3, "asterix"] = "**"
	
	# empty plot (to put horizontal lines before bars), panel.first didn't work with barplot
	barplot(input$d, col="white", border="white", ylim=c(ymin,ymax), width=0.8, space=spaces, xaxt="n", yaxt="n", xlab="")
	
	# horizontal lines
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
		text(x=posi$bars, y=posi$d+posi$se, labels=posi$asterix, pos=3, offset=0, cex=cexnames)
	}
	if(nrow(nega) > 0){
		text(x=nega$bars, y=nega$d-nega$se, labels=nega$asterix, pos=1, cex=cexnames)		
	}
	
	# subtitle
	if (infsites){
		subtitel = paste0(min(input$n_sites)," - ", max(input$n_sites), " informative sites")
		mtext(subtitel, line=-1, adj=0.5, cex=mcex) 
	}

	if (lege[1]){
		# legend Z-scores
		legend("topleft",cex=legcex,ncol=2,legend=c("*", "**", "|Z| > 2","|Z| > 3"),bty="n",title=" weighted block Jackknife")
	}
	if (lege[2]){
		# legend populations
		super = unique(input[,c("super", "color")])
		legend("topright",cex=legcex,ncol=3,legend=super$super,bty="n",fill=colors()[super$color],title="super populations")
	}
}

# replace sample name with official name
get_offi = function(x, offinames){
    offi = offinames[match(x, offinames[,1]), 2]
    offi[is.na(offi)] = x[is.na(offi)]
    return(offi)
}

#### main

if (! interactive()){

	option_list = list(
		make_option(c("-d", "--infile"), type="character", default="out_d", 
			help="d-stats input \n\t\tdefault = %default"),
		make_option(c("-o", "--outpdf"), type="character", default="D_barplot.pdf", 
			help="output pdf file name \n\t\tdefault = %default"),
		make_option(c("-i", "--info"), type="character", 
			help=".csv table with 4 columns: pop, superpop, color, order"),
		make_option(c("-v", "--varpop"), type="integer", default=3,
			help="{1,2,3,4} population that is variable in one plot (x-axis) \n\t\tdefault = %default"),
		make_option(c("-f", "--fixedy"), action="store_true", default=FALSE,
			help="y-axis is fixed for whole input dataset. If not, y-axis is resized for every plot"),
		make_option(c("-a", "--autoswitch"), action="store_true", default=FALSE,
			help="If median D < 0, pop1 and pop2 get switched to always have overall positive values."),
		make_option(c("-n", "--names"), type="character",
			help="optional .csv table with 2 columns: sample name, official name")
	)

	# -i, --info) .csv table with 4 columns like (superpops): 
		# pop			super	color	order
		# CHB			EAS		123		6
		# JPT			EAS		123		6
		# BantuHerero	AFR		91		2
		# BantuKenya	AFR		91		2


	opt_parser = OptionParser(option_list=option_list, description="\nPlot D(pop1, pop2, X, pop4) for every pop1-pop2-pop4 combi as barplot. Position X can be changed. group, order and color bars by info given in .csv file.")
	opt = parse_args(opt_parser)
	print(opt)

	if (! opt$varpop %in% 1:4){
		stop("argument 4 must be 1,2,3 or 4")
	}

	if (is.null(opt$info)){
		stop("Missing info-file -i/--info")
	}

	# read input
	out_d = read.table(opt$infile, as.is=T, header=T)
	superpops = read.csv(opt$info, header=T, as.is=T)

	# pop-combis per page
	fixed_cols = (1:4)[-opt$varpop]
	out_d_pops = out_d[,fixed_cols]
	combis = unique(out_d_pops)
	combis = paste(combis[,1], combis[,2], combis[,3], sep="_")
	out_d_pops = paste(out_d_pops[,1], out_d_pops[,2], out_d_pops[,3], sep="_")

	# switch pop1-pop2 if median D is negative
	if (opt$autoswitch){
		for (i in 1:length(combis)){
			trio_rows = out_d_pops == combis[i]
			if (median(out_d[trio_rows, "d"]) < 0){
				out_d[trio_rows, "d"] = out_d[trio_rows, "d"] * -1
				out_d[trio_rows, "z"] = out_d[trio_rows, "z"] * -1
				out_d[trio_rows,c(1,2)] = out_d[trio_rows,c(2,1)]
			}
		}
	}

	## plot genomewide results
	pdf(opt$outpdf, width=13)
		par(cex.main=0.8, oma=c(3.5,0,0,0), cex.lab=1)
		if (opt$fixedy){
			# compute y-axis range for plot_d_bars (in %)
			out_d$se = out_d$d / out_d$z
			out_d$se[is.na(out_d$se)] = 0  # if D=0 -> Z=0 -> se=NA (0)
			ymin = min(0, min(out_d$d - out_d$se) * 100) * 1.1
			ymax = max(0, max(out_d$d + out_d$se) * 100) * 1.5 # for space above for legends
		} else {
			ymin = NULL
			ymax = NULL
		}
		for (i in 1:length(combis)){
			one_trio = out_d[out_d_pops == combis[i],]
			plot_d_bars(one_trio, superpops, ymin=ymin, ymax=ymax, namesfile=opt$names)
		}
	dev.off()

}
