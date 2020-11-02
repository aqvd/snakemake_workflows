#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

USAGE='USAGE:

	Rsript --vanilla find_unique_peaks.R <MACS2outDir>

> MACS2outDir: directory where MACS2 *.narrowPeak files are. 
			   ¡¡¡ MUST END IN "/" !!!

.narrowPeak files are expected to be named following this pattern:

	{Protein}_{Condition}_{Repeat}_peaks.narrowPeak

{Protein}: the name of the chipped protein
{Condition}: the treatment eg: siC for control, siX for treatment
{Repeat}: the Nth biological repeat of the experiment
	'

## Test arguments are in the right number and format
if (length(args) != 1) {
	stop(USAGE, call.=FALSE)
} else if (length(args) == 1) {
	## Check MACS2 output directory
	PEAKSDIR = args[1]
	print( paste0(">>  ", PEAKSDIR))
	ifelse(grepl("/$", PEAKSDIR),
		print("Correct MACS2outDir name."),
	 	stop('MACS2outDir MUST END IN /',call.=FALSE))

}

print("======= Checks Done ======")


## ================== Load Libraries ==================== ##
library(tidyverse)
library(ChIPseeker)

## ============== Preapare metadata ==================== ##
## Create config table with data for each peak file

configTab <- data.frame(PeakFile=list.files(PEAKSDIR, pattern="narrowPeak",
						full.names=T)) %>%
							mutate_if(is.factor, as.character)



configTab <- configTab %>%
	mutate(Sample=stringr::str_replace(basename(PeakFile),"_peaks.narrowPeak","")) %>%
	mutate(temp=Sample) %>%
	separate(temp, sep = "_", into = c("Prot", "Condition", "Rep", NA)) 
head(configTab)

## ============== Read .narrowPeak ==================== ##
## Create a list with filenames and IDs to read peak files 
samplefiles <- as.list(as.vector(configTab$PeakFile))
names(samplefiles) <- configTab$Sample
samplefiles

peaksList <- lapply(samplefiles,ChIPseeker::readPeakFile)
peaksList[[1]]

##############################################
##		Identify unique peaks per sample	##
##############################################
## Create a list with peaks for each chipped protein in each condition
## list all unique proteins and conditions of the experiment
prots <- as.list(unique(configTab$Prot))
names(prots) <- unique(configTab$Prot)
cond <- as.list(unique(configTab$Condition))
names(cond) <- unique(configTab$Condition)

## Concatenate and merge all peaks per condition for every protein and save
## into grl variable
library(GenomicRanges)

grl <- lapply(prots, function(p){
	a <- lapply(cond, function(c){
		n <- length(peaksList[configTab$Prot==p & configTab$Condition==c])
		cat("P= ", p," C= ", c , " Number Samples= ", n, "\n" )
		
		grL <- GRangesList(peaksList[configTab$Prot==p & configTab$Condition==c])
		merged <- reduce(unlist(grL))
		names(merged) <- c
		return(merged)
		
	})
	#names(a) <- names(a[1])
	return(unlist(a))
})

## Check that the list is named as expected and the content
grl; lapply(grl, function(l) names(l))
## Count the number of peaks for each experiment
npeaks <- lapply(grl, function(l) {
	lapply(l, function(gr) {
		length(gr)
	})
}); unlist(npeaks)

## GRanges operations for each chip
listUniquePeaks <- lapply(grl, function(P) {
	## v_names where the first element is the name of control conditon
	v_conditions <- levels(relevel(factor(names(P)),ref="siC"))
	cont <- v_conditions[1]
	## Prepare lists to save results
	common <- list(v_conditions[-1]); names(common) <- v_conditions[-1]
	contOnly <- list(v_conditions[-1]); names(contOnly) <- v_conditions[-1]
	treatOnly <- list(v_conditions[-1]); names(treatOnly) <- v_conditions[-1]
	## for loop to extract overlaping intervals
	for (treat in v_conditions[-1]) {
		## subsetByOverlaps(x, ranges, invert=c(TRUE,FALSE))
		## subsetByOverlaps returns the subset of x that has an overlap hit
		## with a range in ranges 

		## Common peaks are present in both
		common[[treat]] <- subsetByOverlaps(P[[treat]], P[[cont]])

		## peaks unique for ech condition are the inverse of the overlaps
		contOnly[[treat]] <- subsetByOverlaps(P[[cont]], P[[treat]], invert=TRUE)
		treatOnly[[treat]] <- subsetByOverlaps(P[[treat]], P[[cont]], invert=TRUE)
	}

	return(list(ContOnly = contOnly,
				Common = common, 
				treatOnly = treatOnly)) 

})
## Check that the list is named as expected and the content
listUniquePeaks; lapply(listUniquePeaks, function(P) names(P))
## Count the number of peaks for each experiment
npeaksUnique <- lapply(listUniquePeaks, function(P) {
					lapply(P, function(uniqueIn) {
						lapply(uniqueIn, function(condition) {
							length(condition)
						})
					})
				}); unlist(npeaksUnique)

print("==========NUMBER OF PEAKS:============"); unlist(npeaks); 
print("========NUMBER OF UNIQUE PEAKS:=======");unlist(npeaksUnique)

## Write bed Files with Unique Peaks
l_names <- names(unlist(listUniquePeaks))
condNames <- gsub(x=l_names,pattern='[.]', replacement="_")

## loop for adding unique ID to each peak and write to _uniquePeaks.bed file
for (i in 1:length(l_names) ) {
	## Get peak file
	peaks <- unlist(listUniquePeaks)[[i]] 
	## Uninque ID per peak: "Prot_OnlyIn_Treatment_Number"
	IDname <- gsub(x=condNames[i], pattern='[.]', replacement="_")
	names(peaks) <- paste0(IDname, '_', seq(1:length(peaks)))
	filename <- paste0(PEAKSDIR, "unique/", IDname, "_uniquePeaks.bed")
	rtracklayer::export(peaks, con=filename)
}