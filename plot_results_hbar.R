#!/usr/bin/env Rscript

# Plot D(X, Y, pop3, pop4) for every pop3-pop4 combi as horizontal barplot

library("optparse")

# TODO:
# maybe still use colors for populations as text colors for X and Y

#### helper
  
# compute absolute dynamix x axis maximum 	
dynam_x = function(out_d){
	xmax = max((max(out_d$d + out_d$se)), 0)
	xmin = min(min(out_d$d - out_d$se), 0)
	xmax = max(abs(xmax), abs(xmin))
	xmax = xmax * 1.2
	return(xmax)
}

# replace sample name with official name
get_offi = function(x, offinames){
    offi = offinames[match(x, offinames[,1]), 2]
    offi[is.na(offi)] = x[is.na(offi)]
    return(offi)
}

# plot horizontal barplot for D(pop1, pop2, pop3, pop4) with 3-4 pops fixed combi
plot_d_hbars = function(input, superpops=NULL, xmin=NULL, xmax=NULL, mcex=0.8, legcex=0.7, lege=TRUE, infsites=TRUE, namesfile=NULL){

	# check that pop3 and pop4 are fixed
	pop3 = unique(input$pop3)
	pop4 = unique(input$pop4)
	if (length(pop3) != 1 | length(pop4) != 1){
		print(input)
		stop("One of [pop3, pop4] is not unique.")
	}

	# merge with superpop-info, get order and sort
	if (! is.null(superpops)){
		input$order1 = superpops[match(input$pop1, superpops[,1]), 4]
		input$order2 = superpops[match(input$pop2, superpops[,1]), 4]
		input = input[order(input$order1, input$order2, input$pop1, input$pop2),]
	} else {
		input = input[order(input$pop1, input$pop2),]
	}
	
	# convert to official names
	if (! is.null(namesfile)){
		offinames = read.csv(namesfile, header=T, as.is=T)
		input[,1:4] = apply(input[,1:4], 2, get_offi, offinames)
	}

	# title
	titel = paste0("D(X, Y, ", pop3, ", ", pop4, ")")
	
	# convert to %
	input$se = 100 * (input$d / input$z)
	input$d = 100 * input$d
	
	# replace se=NA (caused by D=0)
	input$se[is.na(input$se)] = 0

	# dynamic x-axis, symmetric
	if (is.null(xmin)) 	{
		xmax = dynam_x(input)
		xmin = -xmax
	}
	span = xmax - xmin
	
	# yaxis, add space on top for legend
	ymin = 0
	ymax = nrow(input) * 1.05
	
	# cex xlab names (0.5-1)
	cexnames = 0.9 - (1/300) * nrow(input)
	cexnames = max(0.5, cexnames)

	# significance asterix
	input$asterix = ""
	input[abs(input$z) > 2, "asterix"] = "*"
	input[abs(input$z) > 3, "asterix"] = "**"
	
	# empty plot (to put horizontal lines before bars), panel.first didn't work with barplot
	barplot(input$d, col="white", border="white", xlim=c(xmin,xmax), ylim=c(ymin,ymax), width=0.8, xaxt="n", yaxt="n", xlab="", horiz=T)

	# horizontal lines
	xaxp = par("xaxp")
	hlines = seq(xaxp[1], xaxp[2], (xaxp[2]-xaxp[1])/xaxp[3])
	abline(v=hlines, col="lightblue", lty=1)

	# barplot
	bars = barplot(input$d, main=titel, xlim=c(xmin,xmax), las=2, width=0.8, xlab="D[%]", cex.axis=0.8, cex.names=cexnames, cex.main=1.2,  add=T, horiz=T, las=1, col="lightblue")

	# box around plot
	box()

	# error bars
	suppressWarnings(arrows(c(input$d+input$se,input$d-input$se), c(bars,bars), c(input$d, input$d), c(bars,bars), angle=90, code=1, length=0.015))

	# pop1,pop2 names
	mtext(side=4, at=bars, text=input$pop1, las=1, line=1, cex=cexnames)
	mtext(side=2, at=bars, text=input$pop2, las=1, line=1, cex=cexnames)

	# X, Y labels
	mtext(side=3, text="Y", adj=-0.1)
	mtext(side=3, text="X", adj=1.1)

	# significance asterix
	input$bars = bars
	offset = 0.02 * span
	ytext = input$d + sign(input$d) * (input$se + offset)
	text(y=input$bars, x=ytext, labels=input$asterix, adj=0.5, cex=cexnames, srt=90)

	# subtitle
	if (infsites){
		subtitel = paste0(min(input$n_sites)," - ", max(input$n_sites), " informative sites")
		mtext(subtitel, line=0.3, adj=0.5, cex=mcex) 
	}

	# legend Z-scores
	if (lege){
		legend("topleft", cex=legcex, ncol=2, legend=c("*", "**", "|Z| > 2","|Z| > 3"), bty="o", bg="white", title=" weighted block Jackknife")
	}
}


#### main
if (sys.nframe() == 0){

	option_list = list(
		make_option(c("-d", "--infile"), type="character", default="out_d", 
			help="d-stats input \n\t\tdefault = %default"),
		make_option(c("-o", "--outpdf"), type="character", default="D_hbarplot.pdf", 
			help="output pdf file name \n\t\tdefault = %default"),
		make_option(c("-i", "--info"), type="character", 
			help="optional .csv table with 4 columns: pop, superpop, color, order. Here for sorting only"),
		make_option(c("-f", "--fixedx"), action="store_true", default=FALSE,
			help="x-axis is fixed for whole input dataset. If not, x-axis is resized for every plot"),
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


	opt_parser = OptionParser(option_list=option_list, description="\nPlot D(X, Y, pop3, pop4) for every pop3-pop4 combi as horizontal barplot")
	opt = parse_args(opt_parser)
	print(opt)

	# read input
	out_d = read.table(opt$infile, as.is=T, header=T)
	if (! is.null(opt$info)){
		superpops = read.csv(opt$info, header=T, as.is=T)
	} else {
		superpops = NULL
	}
	
	# pop-combis per page
	combis = unique(out_d[,c(3,4)])
	combis = paste(combis[,1], combis[,2], sep="_")
	out_d_pops = paste(out_d[,3], out_d[,4], sep="_")

	# switch pop1-pop2 if median D is negative
	if (opt$autoswitch){
		for (i in 1:length(combis)){
			duo_rows = out_d_pops == combis[i]
			if (median(out_d[duo_rows, "d"]) < 0){
				out_d[duo_rows, "d"] = out_d[duo_rows, "d"] * -1
				out_d[duo_rows, "z"] = out_d[duo_rows, "z"] * -1
				out_d[duo_rows,c(1,2)] = out_d[duo_rows,c(2,1)]
			}
		}
	}

	## plot genomewide results
	pdf(opt$outpdf, width=9)
		par(oma=c(0,7,0,9), cex=1, cex.main=1)
		if (opt$fixedx){
			# compute x-axis range, symmetric
			out_d$se = out_d$d / out_d$z
			out_d$se[is.na(out_d$se)] = 0  # if D=0 -> Z=0 -> se=NA (0)
			xmax = dynam_x(out_d)
			xmin = -xmax
		} else {
			xmin = NULL
			xmax = NULL
		}
		for (i in 1:length(combis)){
			one_trio = out_d[out_d_pops == combis[i],]
			plot_d_hbars(one_trio, superpops, xmin=xmin, xmax=xmax, namesfile=opt$names)
		}
	dev.off()

}
