#!/usr/bin/env Rscript

# Plot D(pop1, pop2, X, pop4) per X-freq for every pop1-pop2-pop4 combi as lines or dots
# color by info given in .csv file 


library(optparse)
 
option_list = list(
	make_option(c("-d", "--infile"), type="character", default="out_d_freqbins", 
		help="d-stats input \n\t\tdefault = %default"),
	make_option(c("-o", "--outpdf"), type="character", default="D_freqbins.pdf", 
		help="output pdf file name \n\t\tdefault = %default"),
	make_option(c("-i", "--info"), type="character", 
		help=".csv table with 4 columns: pop, superpop, color, order"),
	make_option(c("-f", "--fixedy"), action="store_true", default=FALSE,
		help="y-axis is fixed for whole input dataset. If not, y-axis is resized for every plot"),
	make_option(c("-s", "--singlepop"), action="store_true", default=FALSE,
		help="Generate one plot per pop3. By default all pop3 for one pop1-pop2-pop4 combi are plotted together."),
	make_option(c("-a", "--autoswitch"), action="store_true", default=FALSE,
		help="If median D < 0, pop1 and pop2 get switched to always have overall positive values."),
	make_option(c("-l", "--lines"), action="store_true", default=FALSE,
		help="Draw lines instead of dots.")
)

# -i, --info) .csv table with 4 columns like (superpops): 
	# pop			super	color	order
	# CHB			EAS		123		6
	# JPT			EAS		123		6
	# BantuHerero	AFR		91		2
	# BantuKenya	AFR		91		2

opt_parser = OptionParser(option_list=option_list, description="\nPlot D(pop1, pop2, X, pop4) per X-freq for every pop1-pop2-pop4 combi as lines or dots. Color by info given in .csv file.")
opt = parse_args(opt_parser)
print(opt)

if (is.null(opt$info)){
	stop("Missing info-file -i/--info")
}



#### helper

# plot D(pop1, pop2, X, pop4) per freqbin
plot_d_freqs = function(input, superpops, ymin=NULL, ymax=NULL){
	
	# titel
	titel = paste0("D(",unique(input$pop1),", ",unique(input$pop2),", X, ",unique(input$pop4),")")
#	subtitel = paste0(min(input$n_sites)," - ", max(input$n_sites), " informative sites")
	ylab="D[%]"
#	xlab="Derived allele frequency in X"
	xlab=""

	# merge with superpop-info and define spaces between bars
	plot_pops = unique(input$pop3)
	plot_cols = colors()[superpops[match(plot_pops, superpops[,1]), 3]]

	# empty plot
	# TODO: make xlim parameter for zoom-ins
	plot(input$bin, input$d, type="n", xlim=c(0,1.07), ylim=c(ymin,ymax), ylab=ylab, xlab=xlab)  
	# add grid lines
	abline(h=0, col="lightblue")
	for (k in seq(0.01, 0.05, 0.01)){
		abline(v=k, col="lightgray")
	}
	# titel
	mtext(titel, line=2, adj=0.5, cex=par()$cex.main, font=2)
#	mtext(subtitel, line=-2, adj=0.5, cex=par()$cex.lab) 
	# legend
	legend("topright", legend=plot_pops, lty=1, col=plot_cols, bty="n", cex=1.2, lwd=3)
	# add points/lines for each pop3
	for (j in 1:length(plot_pops)){
		one3 = input[input$pop3 == plot_pops[j],]
		if(opt$lines){
			lines(one3$bin, one3$d, col=plot_cols[j])
			points(one3$bin, one3$d, col=plot_cols[j], cex=0.5)
		} else {
			points(one3$bin, one3$d, col=plot_cols[j], cex=0.7)
		}
	}
}

# plot range of # of informative sites for each freqbin
plot_n_sites = function(input){
	xlab="Derived allele frequency in X"
	ylab="# ABBA/BABA sites"
	subtitel = paste0(min(input$n_sites)," - ", max(input$n_sites), " informative sites")
	
	mean_sites = tapply(input$n_sites, input$bin, mean) 
	min_sites = tapply(input$n_sites, input$bin, min) 
	max_sites = tapply(input$n_sites, input$bin, max)
	x = as.numeric(names(min_sites)) # TODO: make xlim parameter for zoom-ins
	plot(x, mean_sites, type="h", lwd=4, col="steelblue", ylim=c(0, max(max_sites)), xlim=c(0,1.05)
		, xlab="", ylab=ylab)
	segments(x0=x, y0=min_sites, y1=max_sites)
	
	mtext(xlab, side=1, line=3.5, adj=0.5, cex=par()$cex.lab-0.2)
	mtext(subtitel, line=-2, adj=0.5, cex=par()$cex.sub-0.2)
	
}


#### main

# read input
d_freqs = read.table(opt$infile, as.is=T, header=T)
superpops = read.csv(opt$info, header=T, as.is=T)

# convert D to percent
d_freqs$d = d_freqs$d * 100

# pop-combi per page
if (opt$singlepop){
	d_freqs$paste_pop = paste(d_freqs[,1], d_freqs[,2], d_freqs[,3], d_freqs[,4], sep="_")
} else {
	d_freqs$paste_pop = paste(d_freqs[,1], d_freqs[,2], d_freqs[,4], sep="_")
}
combis = unique(d_freqs$paste_pop)


# switch pop1-pop2 if median D is negative
if (opt$autoswitch){
	for (i in 1:length(combis)){
		page_rows = d_freqs$paste_pop == combis[i]
		if (median(d_freqs[page_rows, "d"]) < 0){
			d_freqs[page_rows, "d"] = d_freqs[page_rows, "d"] * -1
			d_freqs[page_rows,c(1,2)] = d_freqs[page_rows,c(2,1)]
		}
	}
}	


## plot page per pop-combi
pdf(opt$outpdf, width=11, height=8)
	par(cex.main=1.2, cex.lab=1.2, oma=c(5,1.5,4,1), mar=c(1,4,2,2))
	layout(matrix(c(1,1,2),3,1))
	if (opt$fixedy){
		# compute y-axis range for plot_d_freqs (in %)
		ymin = min(0, min(d_freqs$d)) * 1.1
		ymax = max(0, max(d_freqs$d)) * 1.1 # for space above for legends
	} else {
		ymin = NULL
		ymax = NULL
	}
	for (pp in combis){
		page_data = d_freqs[d_freqs$paste_pop == pp,]
		plot_d_freqs(page_data, superpops, ymin=ymin, ymax=ymax)
		plot_n_sites(page_data)
	}
dev.off()

