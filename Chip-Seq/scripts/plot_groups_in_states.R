#!/usr/bin/env Rscript

if (!require("pacman")) install.packages("pacman", repos="https://cran.rediris.es")
pacman:::p_load(tidyverse, optparse)

option_list = list(
 	optparse::make_option(c("-d", "--dir"), type="character",
 				default=NULL, metavar="character",
				help="directory with dorted .bed divided in groups(--dir)"),

 	optparse::make_option(c("-r", "--regionsRegex"), type="character",
				default="_log2ratio.bed$", metavar="character",
				help="Regex to identify sorted .bed files"),

 	optparse::make_option(c("--protCondRegex"), type="character",
				default="(.+)_unique.+_(.+)_vs_.+[.]bed_(.+)", metavar="character",
				help="Regex to extract Protein Condition and State from .bed files"),

 	optparse::make_option(c("-O", "--outDir"), type="character",
    			default=NULL, metavar="character",
				help="output directory created inside --dir" ),

 	optparse::make_option(c("--width"), type="numeric",
    			default=7, metavar="integer",
				help="width of plot in inches" ),

 	optparse::make_option(c("--height"), type="numeric",
    			default=7, metavar="integer",
				help="height of plot in inches" )
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dir)){
	print_help(opt_parser)
	stop("Specify directory with dorted .bed divided in groups(--dir).", call.=FALSE)
} else if (is.null(opt$regionsRegex)) {
	print_help(opt_parser)
	stop("Specify regex to identify -bedFiles (--regionsRegex).", call.=FALSE)
} else if (is.null(opt$outDir)) {
	print_help(opt_parser)
	stop("Specify output directory created inside --dir (--outDir).", call.=FALSE)
} else if (is.null(opt$protCondRegex)) {
	print_help(opt_parser)
	stop("Specify regex to extract prot cond from .bed files (--protCondRegex).", call.=FALSE)
}

## Get values from options
wd = opt$dir
regionsRegex = opt$regionsRegex
outdir = opt$outDir
width = opt$width
height = opt$height
protCondRegex = opt$protCondRegex

## Tests! remove
# wd <- "/Users/aqo/Desktop/cluster/mouseHepatocytes_NIPBLko/log2Ratio_NIPBLko_vs_WT/10groups"
# regionsRegex <- "_CTCF$|bed_promoters$"
# outdir <- file.path(wd, "plots_groups_in_states")

## Create wd
dir.create(outdir, showWarnings = FALSE)
setwd(outdir)

## ================================================================== ##
## 								MAIN 								  ##
## ================================================================== ##

## List files
groupFiles <- list.files(wd, pattern = regionsRegex)
res <- data.frame()

for (i in groupFiles){
	prot_cond_state <- stringr::str_match(i, protCondRegex) %>%
		data.frame() %>%
		dplyr::mutate(prot= gsub("_.+", "", X2), cond= X3, 
			state= gsub(".+[_]?.+_(.+)", "\\1", X4))
	print(prot_cond_state)

	bed <- read.table(file.path(wd, i), sep="\t", header=FALSE)
	toRbind <- data.frame(group= bed$V4, prot= prot_cond_state[1,"prot"],
		cond= prot_cond_state[1,"cond"], state=prot_cond_state[1,"state"])

	res <- rbind(res, toRbind)
}
res <- as_tibble(res)
res

plotPath <- file.path(wd, outdir, "barplot_CTCF_PROM_2columns.png")
print(plotPath)
ggplot(res) + 
	geom_bar(aes(as.factor(group),fill=group),alpha=0.8) + 
	scale_fill_gradient(low="red", high="blue") + 
	facet_wrap(~prot + cond + state, ncol=2, scales="free") + 
	xlab("Group") + 
	ggsave(plotPath, height = height, width= width)






