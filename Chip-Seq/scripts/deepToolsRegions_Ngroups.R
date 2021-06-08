#!/usr/bin/env Rscript

if (!require("pacman")) install.packages("pacman", repos="https://cran.rediris.es")
pacman:::p_load(tidyverse, optparse)

option_list = list(
 	optparse::make_option(c("-d", "--dir"), type="character",
 				default=NULL, metavar="character",
				help="Directory with log2FC sorted regions"),

 	optparse::make_option(c("-r", "--regionsRegex"), type="character",
				default="_log2ratio.bed$", metavar="character",
				help="Regex to identify sorted .bed files"),

 	optparse::make_option(c("-n", "--nGroups"), type="integer",
				default=10, metavar="character",
				help="Number of groups to divide regions in .bed files"),

    optparse::make_option(c("-O", "--outDir"), type="character",
    			default=NULL, metavar="character",
				help="output directory" )

); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dir)){
	print_help(opt_parser)
	stop("Specify directory (--dir).n", call.=FALSE)
} else if (is.null(opt$regionsRegex)) {
	print_help(opt_parser)
	stop("Specify regex to identify -bedFiles (--regionsRegex).n", call.=FALSE)
} else if (! is.integer(opt$nGroups)) {
	print_help(opt_parser)
	stop("Integer many groups to divide .bed (--nGroups).n", call.=FALSE)
} else if (is.null(opt$outDir)) {
	print_help(opt_parser)
	stop("Specify output directory(--outDir).n", call.=FALSE)
}

## Get values from options
wd = opt$dir
regionsRegex = opt$regionsRegex
NGROUPS = opt$nGroups
outdir = opt$outDir

# wd <- "/Users/aqo/Desktop/cluster/HAP1_cells_MAU2ko_WAPLko/log2Ratio_MAU2ko_and_WAPLko_vs_WT"
# NGROUPS <- 10
# regionsRegex <- "_log2ratio.bed$"
# outdir <- file.path(wd, "10groups_header")
## Create wd
dir.create(outdir, showWarnings = FALSE)
setwd(outdir)

getwd()

## ================================================================== ##
## 								MAIN 								  ##
## ================================================================== ##

## List files
log2files <- list.files(wd, pattern = regionsRegex)

## read files into tables
bed_l <- lapply(log2files, function(file){
	bed <- read.table(file.path(wd,file), sep="\t", header = TRUE, 
		comment.char="") %>%
		as_tibble() %>%
		`names<-`(gsub("X.","", colnames(.)))

	t = ceiling( ((nrow(bed)) / NGROUPS))

	groups <- rep(1:NGROUPS,times=1, each=t)[1:nrow(bed)]
	length(groups) == nrow(bed)

	bed <- bed %>% dplyr::mutate(group=as.factor(groups)) %>%
		dplyr::select(chrom, start, end, group) %>%
		dplyr::rename(`#chrom`="chrom") # valid header for bedtools needs "#"
	
	return(bed)
}); names(bed_l) <- make.names(log2files)

## Summarize number of peaks per group per state
lapply(bed_l, function(bed){
	bed %>% group_by(group) %>%
		dplyr::summarize(n=n()) 
})

## Write bed files 
sapply(names(bed_l), function(i){
	out_filename <- paste0(i, "_", NGROUPS, "_groups.bed")
	write.table(bed_l[[i]], file=out_filename, quote = FALSE, 
		row.names = FALSE, col.names = TRUE, sep="\t")
})












