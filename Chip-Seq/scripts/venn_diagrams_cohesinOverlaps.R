library(tidyverse)
library(eulerr)

## ================================================================== ##
## This first commented part was done by hand using bedtools, but i had some
## more number of peaks after intersections than the original files. 
## So I don't know if trust that. The next part using mergePeaks give similar
## results with respect to the number of peaks, SMC1, intersected. So keep using
## that.
## ================================================================== ##


# wd="/Users/aqo/Desktop/MCF10A/siNipbl/Peaks/merged_replicates_callpeak/intersected_with_SMC1/SA1only_common_SA2only"
# setwd(wd)

# files <- list.files(wd, pattern = "toVenn")
# toEuler <- sapply(files, function(f){
# 	read.table(f)
# })

# names(toEuler) <- names(toEuler) %>% stringr::str_replace(.,"_toVenn.+", "")
# lapply(toEuler, length)

# ix <- grep("SMC1", names(toEuler))
# cols <- c("lightblue2", "skyblue3", "sienna1", "red4", "lightgreen", "limegreen")
# filename <- paste0(wd, "/euler_SA1-SA2_intersected_SMC1.pdf")
# p1 <- plot(euler(toEuler[-ix]),
# 	quantities = list(font=3, fontsize = 8),
# 	fills= list(alpha=0.8,fill=cols[-ix]),
# 	labels=list(TRUE, labels=gsub("[_]", " ", names(toEuler)),
# 		fontsize=6),
# 	edges=FALSE,
# 	legend=list(TRUE,
# 		fontsize=8),
# 	main=list(label="STAG peaks overlaping SMC1 peaks",
# 		font=3,  fontsize = 10)) 

# pdf(filename, width=6.5, height=3); print(p1); dev.off()

## ================================================================== ##
## 						mergePeaks process output					  ##
## ================================================================== ##
# wd="/Users/aqo/Desktop/MCF10A/siNipbl/Peaks/merged_replicates_callpeak/SAonly-common_allPeaks"
wd="/Users/aqo/Desktop/cluster/MCF10A/siNipbl/Peaks/merged_replicates_callpeak/SAonly-common_allPeaks/Cohesin_intersected_CTCF"

## Read Homer merge peaks venn file
files <- list.files(wd, pattern = "mergePeaks_venn.txt$", full.names = T)

toEuler <- lapply(files, function(f){
	## Remove full path form column names
	venn_table_df<-read.table(f,header = TRUE,sep = "\t",
			stringsAsFactors = FALSE,check.names = FALSE) %>%
		`colnames<-`(gsub("/.+/","", colnames(.))) %>%
		dplyr::mutate(Name=stringr::str_split(Name,"\\|") %>% 
		purrr::map(., ~stringr::str_replace(.,"/.+/","") %>% 
			gsub("_merged_peaks.narrowPeak", "", .)) %>%
		purrr::map(., ~paste(., collapse = "&")))

	toEuler<-venn_table_df[['Total']]
	names(toEuler) <- venn_table_df[['Name']]

	return(toEuler)
}); names(toEuler) <- gsub("/.+/","",files)


for ( venn_filename in names(toEuler) ) {
	overlaps <- toEuler[[venn_filename]]
	
	## Fit euler and write residuals (some regions can be missed)
	fit <- euler(overlaps, shape="ellipse")

	resid_filename <- file.path(wd, paste0(venn_filename, "_residuals.tsv"))
	options(scipen=999)
		fit$residuals %>% data.frame() %>% 
			write.table(file=resid_filename, sep="\t", col.names=F, quote=FALSE)
	options(scipen=0)

	p1 <- plot(fit,
		quantities = list(font=4, fontsize = 5),
		fills= list(alpha=0.8, 
			# fill=c("limegreen",
			# 	"red4",
			# 	"skyblue3",
			# 	"sienna1",
			# 	"lightgreen",
			# 	"sienna3",
			# 	"plum2")
			# 	),
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
	filename <- paste0(wd, "/euler_", f, ".pdf")

	# save the plot into a PDF
	pdf(filename, width=6.5, height=3); print(p1); dev.off()
}

## ============= Get coordenates of those peaks ========= ##
mergePeaks_files <- list.files(wd, pattern = "mergePeaks.txt$") 

mergePeaks_l <- plyr::llply(mergePeaks_files, .progress='text', function(f){
	read.table(file.path(wd, f), sep="\t", header=T,
		stringsAsFactors=F) %>% as.tibble()
}); names(mergePeaks_l) <- mergePeaks_files

plyr::llply(mergePeaks_l, function(x){
	nested <- x %>% dplyr::group_by(Parent.files) %>% tidyr::nest() %>%
		mutate(Parent.files=stringr::str_split(Parent.files,"\\|") %>% 
			purrr::map(., ~stringr::str_replace(.,"/.+/","") %>% 
				gsub("_merged_peaks.narrowPeak", "", .)) %>%
			purrr::map_chr(., ~paste(., collapse = "@"))) %>%
		# dplyr::filter(grepl("@SMC1|SMC1@", Parent.files)) %>% To keep SMC1
		dplyr::mutate(bed=purrr::map(data, ~ dplyr::select(., chr,start,end)))

	apply(nested, 1, function(i){
		write.table(i["bed"], file=paste0(wd,"/", i["Parent.files"], ".bed"),
				col.names=FALSE, sep="\t", quote=FALSE, row.names = FALSE)
	})

})

## =========== CHECK THAt WE DO NOT LOOSE ANY PEAK USING MERGE PEAKS ======== ##
## Read siC peaks and create nested data frame with rows containing peaks
## in venn intersection 
nested <- mergePeaks_l[[1]] plyr::group_by(Parent.files) %>% tidyr::nest() %>%
		mutate(Parent.files=stringr::str_split(Parent.files,"\\|") %>% 
			purrr::map(., ~stringr::str_replace(.,"/.+/","") %>% 
				gsub("_merged_peaks.narrowPeak", "", .)) %>%
			purrr::map_chr(., ~paste(., collapse = "@"))) %>%
		# dplyr::filter(grepl("@SMC1|SMC1@", Parent.files)) %>% To keep SMC1
		dplyr::mutate(bed=purrr::map(data, ~ dplyr::select(., chr,start,end)))

## Columns starting by X.Users (path to peak file) contain the IDs of the called
## peaks in each of the merged peak. So Extract those IDs, find unique ones and sum
nested %>% dplyr::mutate(IDs=purrr::map(
	data, ~ 
		(.) %>% tidyr::pivot_longer(starts_with("X.users"), names_to="Peak_file",
			names_prefix="X.+replicates_callpeak[.]", values_to="ID") %>%
			tidyr::separate_rows(ID, sep=",") %>% pull(ID)
		)
	) %>%
	dplyr::mutate(unique_peaks=purrr::map(IDs, ~unique(.))) %>%
	pull(unique_peaks) %>% unlist(., use.names = F) %>% unique %>% length

	