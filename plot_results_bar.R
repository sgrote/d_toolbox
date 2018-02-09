#!/usr/bin/env Rscript

# Plot D(pop1, pop2, X, pop4) for every pop1-pop2-pop4 combi as barplot
# group, order and color bars by info given in .csv file 


library("optparse")
 
option_list = list(
	make_option(c("-d", "--infile"), type="character", default="out_d", 
		help="d-stats input"),
	make_option(c("-o", "--outpdf"), type="character", default="D_barplot.pdf", 
		help="output pdf file name"),
	make_option(c("-i", "--info"), type="character", 
		help=".csv table with 4 columns: pop, superpop, color, order"),
	make_option(c("-v", "--varpop"), type="integer", default=3,
		help="{1,2,3,4} population that is variable in one plot (x-axis)"),
	make_option(c("-f", "--fixedy"), action="store_true", default=FALSE,
		help="y-axis is fixed for whole input dataset. If not, y-axis is resized for every plot")
)

# -i, --info) .csv table with 4 columns like (superpops): 
	# pop			super	color	order
	# CHB			EAS		123		6
	# JPT			EAS		123		6
	# BantuHerero	AFR		91		2
	# BantuKenya	AFR		91		2

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print(opt)

if (! opt$varpop %in% 1:4){
	stop("argument 4 must be 1,2,3 or 4")
}


# source d_functions.R  (plot_d_bars())
argv = commandArgs(trailingOnly = FALSE)
script_dir = dirname(substring(argv[grep("--file=", argv)], 8))
source(paste0(script_dir, "/d_functions.R"))


##############

# read input
out_d = read.table(opt$infile, as.is=T, header=T)
superpops = read.csv(opt$info, header=T, as.is=T)

# pop-combis per page
fixed_cols = (1:4)[-opt$varpop]
out_d_pops = out_d[,fixed_cols]
combis = unique(out_d_pops)
combis = paste(combis[,1], combis[,2], combis[,3], sep="_")
out_d_pops = paste(out_d_pops[,1], out_d_pops[,2], out_d_pops[,3], sep="_")

## plot genomewide results
pdf(opt$outpdf, width=13)
	par(cex.main=0.8, oma=c(2.5,0,0,0), cex.lab=1)
	if (opt$fixedy){
		# compute y-axis range for plot_d_bars (in %)
		out_d$se = out_d$d / out_d$z
		out_d$se[is.na(out_d$se)] = 0  # if D=0 -> Z=0 -> se=NA (0)
		ymin = min(0, min(out_d$d - out_d$se) * 100) * 1.1
		ymax = max(0, max(out_d$d + out_d$se) * 100) * 1.4 # for space above for legends
		sw = FALSE  # don't switch pop1 and pop2 for overall negative D (fixed ylim not appropriate for that)
	} else {
		ymin = NULL
		yamx = NULL
		sw = TRUE
	}
	for (i in 1:length(combis)){
		one_trio = out_d[out_d_pops == combis[i],]
		plot_d_bars(one_trio, superpops, ymin=ymin, ymax=ymax, auto_switch=sw)
	}
dev.off()

