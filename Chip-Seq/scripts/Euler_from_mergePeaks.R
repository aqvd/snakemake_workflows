#!/usr/bin/env Rscript

## First specify the packages of interest

if (!require("pacman")) install.packages("pacman", repos="https://cran.rediris.es")
pacman:::p_load(tidyverse, eulerr, optparse)

 
option_list = list(
 	optparse::make_option(c("-v", "--vennFileRegex"), type="character",
 				default="mergePeaks_venn.txt$", metavar="character",
				help="HOMER merge peaks file generated wit --venn"),

 	optparse::make_option(c("-f", "--mergedPeaksRegex"), type="character",
				default="mergePeaks.txt$", metavar="character",
				help="file that has HOMER mergePeaks output"),

	optparse::make_option(c("-d", "--mergePeaksDir"), type="character",
				default=NULL, metavar="character",
				help="directory where HOMER mergePeaks files are"),

    optparse::make_option(c("-O", "--outDir"), type="character",
    			default=NULL, metavar="character",
				help="output directory" )

); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$vennFileRegex)){
	print_help(opt_parser)
	stop("At least one argument must be supplied (vennFileRegex).n", call.=FALSE)
} else if (is.null(opt$mergedPeaksRegex)) {
	print_help(opt_parser)
	stop("At least one argument must be supplied (mergedPeaksRegex).n", call.=FALSE)
} else if (is.null(opt$outDir)) {
	print_help(opt_parser)
	stop("At least one argument must be supplied (outDir).n", call.=FALSE)
} else if (is.null(opt$mergePeaksDir)) {
	print_help(opt_parser)
	stop("At least one argument must be supplied (outDir).n", call.=FALSE)
}

vennFileRegex = opt$vennFileRegex
mergedPeaksFile = opt$mergedPeaksRegex
outdir = opt$outDir
mergePeaks_dir = opt$mergePeaksDir

## ================================================================== ##
## 						mergePeaks process output					  ##
## ================================================================== ##

files <- list.files(mergePeaks_dir, pattern = "mergePeaks_venn.txt$",
	full.names = T)

toEuler <- lapply(files, function(f){
	## Remove full path form column names
	venn_table_df<-read.table(f,header = TRUE,sep = "\t",
			stringsAsFactors = FALSE,check.names = FALSE) %>%
		`colnames<-`(gsub("/.+/","", colnames(.))) %>%
		dplyr::mutate(Name=stringr::str_split(Name,"\\|") %>% 
		purrr::map(., ~stringr::str_replace(.,"/.+/","") %>% 
			gsub("_peaks.narrowPeak", "", .)) %>%
		purrr::map(., ~paste(., collapse = "&")))

	toEuler<-venn_table_df[['Total']]
	names(toEuler) <- venn_table_df[['Name']]

	return(toEuler)
}); names(toEuler) <- gsub("/.+/","",files)
cols <- c("lightblue2", "skyblue3", "sienna1", "red4","lightgreen", "limegreen")

for ( venn_filename in names(toEuler) ) {
	overlaps <- toEuler[[venn_filename]]
	
	p1 <- plot(euler(overlaps),
	quantities = list(font=4, fontsize = 5),
	fills= list(alpha=0.8, 
				fill=c(
				"sienna1",
				"limegreen",
				"red4",
				"skyblue3",
				"lightgreen",
				"sienna3",
				"plum2")
				),
	labels=list(TRUE, labels=names(overlaps)[-grep("&",names(overlaps))],
		fontsize=8),
	edges=FALSE,
	legend=list(TRUE,
		fontsize=6),
	main=list(label='All "Cohesin" peaks',
		font=3,  fontsize = 10)) 

	f <- gsub("/.+/(.+)_mergePeaks_venn.txt", "\\1", venn_filename)
	filename <- paste0(outdir, "/euler_", f, ".pdf")

	# save the plot into a PDF
	pdf(filename, width=6.5, height=3); print(p1); dev.off()
}

## ============= Get coordenates of those peaks ========= ##
mergePeaks_files <- list.files(mergePeaks_dir, pattern = "mergePeaks.txt$") 

mergePeaks_l <- plyr::llply(mergePeaks_files, .progress='text', function(f){
	read.table(file.path(mergePeaks_dir, f), sep="\t", header=T,
		stringsAsFactors=F) %>% as.tibble()
}); names(mergePeaks_l) <- mergePeaks_files

plyr::llply(mergePeaks_l, function(x){
	nested <- x %>% dplyr::group_by(Parent.files) %>% tidyr::nest() %>%
		mutate(Parent.files=stringr::str_split(Parent.files,"\\|") %>% 
			purrr::map(., ~stringr::str_replace(.,"/.+/","") %>% 
				gsub("_peaks.narrowPeak", "", .)) %>%
			purrr::map_chr(., ~paste(., collapse = "@"))) %>%
		# dplyr::filter(grepl("@SMC1|SMC1@", Parent.files)) %>% To keep SMC1
		dplyr::mutate(bed=purrr::map(data, ~ dplyr::select(., chr,start,end)))

	apply(nested, 1, function(i){
		write.table(i["bed"], file=paste0(outdir,"/", i["Parent.files"], ".bed"),
				col.names=FALSE, sep="\t", quote=FALSE, row.names = FALSE)
	})

})

