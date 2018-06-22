#!/usr/bin/env Rscript

# plot ABBA and BABA counts for D(pop1, pop2, pop3, pop4) per ANCESTRAL-DERIVED allele pair
# takes output from 'abba_baba_allele_counts.R' as input

library(optparse)
library(plyr)

option_list = list(
	make_option(c("-d", "--infile"), type="character", default="abba_baba_per_base",
		help="d-stats input \n\t\tdefault = %default"),
	make_option(c("-o", "--outpdf"), type="character", default="abba_baba_per_base.pdf",
		help="output file name \n\t\tdefault = %default")
)

opt_parser = OptionParser(option_list=option_list, description="\nplot ABBA and BABA counts for D(pop1, pop2, pop3, pop4) per ANCESTRAL-DERIVED allele pair.\n
takes output from 'abba_baba_allele_counts.R' as input")
opt = parse_args(opt_parser)
print(opt)


## helper
# plot one pop-combi
plotty = function(county){
	par(cex.axis=1.2)
	county$d = (county$baba - county$abba) / (county$baba + county$abba) * 100
	county$d[is.na(county$d)] = 0
	d_total = sum(county$baba - county$abba) / sum(county$baba + county$abba) * 100
	
	titel = paste0("D(", paste(unique(county[,1:4]), collapse=", "),")")
	subtitel = paste0("D = ", round(d_total, digits=2), "%")
	
	mat = as.matrix(county[,c("abba","baba")])
	row.names(mat) = paste(county$A_allele, county$B_allele, sep="_")
	
	ylim1 = c(0, max(c(county$abba, county$baba)) * 1.3)
	d_range = abs(max(county$d) - min(county$d))
	ylim2 = c(min(c(0,county$d)) - 0.1 * d_range, max(c(0,county$d)) + 0.3 * d_range)
	
	layout(matrix(c(1,1,2),1,3))
	barplot(mat, beside=T, col=coly, ylim=ylim1,
		names.arg=c("ABBA","BABA"), legend=T,
	args.legend=list(title="A_B alleles", bty="n", ncol=2, x="topleft", cex=1.3))
	barplot(county[,"d"], col=coly, names.arg="D[%]", ylim=ylim2)
	mtext(titel, 3, -3, outer=T)
	mtext(subtitel, 3,-5, outer=T)
}

## main
counts = read.table(opt$infile, as.is=T, header=T)
coly = colors()[c(43,46,48,51,33,36,148,151,103,106,134,137)]
coly = coly[1:nrow(unique(counts[,c("A_allele", "B_allele")]))]

pdf(opt$outpdf, width=9)
	nothing = ddply(counts, .(pop1,pop2,pop3,pop4), plotty)
dev.off()
