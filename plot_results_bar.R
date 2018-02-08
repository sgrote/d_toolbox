
# Plot D(pop1, pop2, X, pop4) for every pop1-pop2-pop4 combi as barplot
# group, order and color bars by info given in .csv file 

### command line arguments:
# 1) infile (out_d)
# 2) outfile.pdf
# 3) .csv table with 4 columns like (superpops): 
	# pop			super	color	order
	# CHB			EAS		123		6
	# JPT			EAS		123		6
	# BantuHerero	AFR		91		2
	# BantuKenya	AFR		91		2
# 4) T/TRUE for fixed y axis (min and max across whole infile)

args = commandArgs(trailingOnly=TRUE)
out_d_file = args[1]
out_pdf = args[2]
superpops_file = args[3]
if (length(args)==4){
	fixed_y = args[4]
} else {
	fixed_y = FALSE
}
######################################


source("/mnt/expressions/steffi/D/d_toolbox/d_functions.R") # plot_d_bars

out_d = read.table(out_d_file, as.is=T, header=T)
superpops = read.csv(superpops_file, header=T, as.is=T)


## plot genomewide results
pdf(out_pdf, width=13)
	par(cex.main=0.8, oma=c(2.5,0,0,0), cex.lab=1)
	combis = unique(out_d[,c("pop1", "pop2", "pop4")])
	if (fixed_y){
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
	for (i in 1:nrow(combis)){
		one124 = out_d[out_d$pop1==combis[i, "pop1"] & out_d$pop2==combis[i,"pop2"] & out_d$pop4==combis[i,"pop4",],]
		plot_d_bars(one124, superpops, ymin=ymin, ymax=ymax, auto_switch=sw)
	}
dev.off()

