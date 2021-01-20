import os
import  pandas as pd
import numpy as np
import glob
import re

#UNIQUEDIR="/storage/scratch01/users/aquevedo/snakemakeTest/res/macs/unique/"
# UNIQUEDIR="/Users/aqo/Desktop/siNipblChipPeaks/unique/"
UNIQUEDIR="/Users/aqo/Desktop/MCF10A/siNipbl/Peaks/merged_replicates_callpeak/unique/"

## Data frame with .narrowPeak files
naPeak = glob.glob(UNIQUEDIR + "*uniquePeaks.bed")
meta=pd.DataFrame(naPeak, columns=["uniqPeakFile"])

## Extract basename
meta["Basename"]=[re.sub(pattern='^.+/', repl="",string=f) for f in meta.uniqPeakFile]

## EXtract Prot Condition and Rep from Basename
meta[["Protein","Condition","uniqueIn","tmp"]]=meta["Basename"].str.split("_",expand=True)
meta=meta.drop(columns=['tmp'])

meta['Sample']=meta.Protein +"_"+meta.uniqueIn+"_"+ meta.Treatment

## add Background col to know which file use as bg in homer
# meta['Background']= [ meta.uniqPeakFile[(meta.Protein == P) & 
# 					(meta.uniqueIn == 'Common')].values[0] \
#                      if (UniqIn != 'Common') else "" \
#                      for P,UniqIn in zip(meta.Protein, meta.uniqueIn) ]
## >>>> NEW! added after merging replicates in siNipbl chip >>>
# meta["Common"] = ["yes" if 'Common' in UniqIn 
# 					else ""
# 					for UniqIn in meta.uniqueIn]
#
# meta['Background'] = [meta.uniqPeakFile[
# 								(meta.Protein == Prot) & 
# 								(meta.Common == 'yes') &
# 								(meta.Treatment == Treat) ].values[0] \
#                      if (Comm != "yes") else "" \
#                      for Prot,Treat,Comm in zip(meta.Protein,meta.Treatment, meta.Common)]
## <<<< NEW! added after merging replicates in siNipbl chip <<<<<
print(meta.uniqPeakFile.values)
print("===========================")
print(meta.Background.values)

############################################################################
############					RULES							############	
############################################################################

rule all:
	input:
		# expand(UNIQUEDIR + '{Sample}_unique.finished.txt', 
		# 	Sample=meta.Sample[meta.uniqueIn != 'Common'].unique()),
		# expand(UNIQUEDIR + '{Sample}_common.finished.txt', 
		# 	Sample=meta.Sample[meta.uniqueIn == 'Common'].unique()),
		# expand(UNIQUEDIR + 'annot/{Sample}_uniq.txt',
		# 	Sample=meta.Sample[meta.uniqueIn != 'Common'].unique()),
		# expand(UNIQUEDIR + 'annot/{Sample}_common.txt',
		# 	Sample=meta.Sample[meta.uniqueIn == 'Common'].unique()),
		## >>>>>>>>>>
		expand(UNIQUEDIR + '{Sample}_unique.finished.txt', 
			Sample=meta.Sample[meta.Common != 'yes'].unique()),
		expand(UNIQUEDIR + '{Sample}_common.finished.txt', 
		 	Sample=meta.Sample[meta.Common == 'yes'].unique()),
		## <<<<<<<<<<

rule homer_motifs_unique:
	input:
		# peaks= lambda wildcards: expand("{UniqBED}",
		# 	UniqBED=meta.uniqPeakFile[
		# 	(meta.uniqueIn != 'Common') & (meta.Sample == wildcards.Sample)]),
		# background=lambda wildcards: expand("{bgBED}",
		# 	bgBED=meta.Background[(meta.uniqueIn != 'Common') &
		# 						(meta.Sample == wildcards.Sample)])
		peaks= lambda wildcards: expand("{UniqBED}",
			UniqBED=meta.uniqPeakFile[
			(meta.Common != 'yes') & (meta.Sample == wildcards.Sample)]),
		background=lambda wildcards: expand("{bgBED}",
			bgBED=meta.Background[
			(meta.Common != 'yes') & (meta.Sample == wildcards.Sample)])
	output:
		UNIQUEDIR + '{Sample}_unique.finished.txt'
	params:
		outdir=UNIQUEDIR + '{Sample}_homerMotifs'
	log:
		UNIQUEDIR+'log/{Sample}_unique.log'
	threads: 1
	shell:
		'findMotifsGenome.pl {input.peaks} hg19 {params.outdir} \
		-bg {input.background} -size given &> {log} && touch {output}'

rule homer_motifs_common:
	input:
		peaks= lambda wildcards: expand("{CommonBED}",
			CommonBED=meta.uniqPeakFile[
			(meta.uniqueIn == 'Common') & (meta.Sample == wildcards.Sample)]),
	output:
		UNIQUEDIR + '{Sample}_common.finished.txt'
	params:
		outdir=UNIQUEDIR + '{Sample}_homerMotifs'
	log:
		UNIQUEDIR+'log/{Sample}_common.log'
	threads: 1
	shell:
		'findMotifsGenome.pl {input.peaks} hg19 {params.outdir} \
		-size given &> {log} && touch {output}'

rule homer_annotatePeaks_unique:
	input:
		peaks= lambda wildcards: expand("{UniqBED}",
			UniqBED=meta.uniqPeakFile[
			(meta.uniqueIn != 'Common') & (meta.Sample == wildcards.Sample)])
	output:
		anno=UNIQUEDIR + 'annot/{Sample}_uniq.txt',
		stats=UNIQUEDIR + 'annot/{Sample}_uniq_stats.txt'
	threads: 2
	params:
		genome="hg19",
	log:
		UNIQUEDIR+'log/{Sample}_annoUniq.log'
	shell:
		'annotatePeaks.pl {input.peaks} {params.genome} -size given\
		-annStats {output.stats} -cpu {threads} 1> {output.anno} 2> {log}'

rule homer_annotatePeaks_common:
	input:
		peaks= lambda wildcards: expand("{UniqBED}",
			UniqBED=meta.uniqPeakFile[
			(meta.uniqueIn == 'Common') & (meta.Sample == wildcards.Sample)])
	output:
		anno=UNIQUEDIR + 'annot/{Sample}_common.txt',
		stats=UNIQUEDIR + 'annot/{Sample}_common_stats.txt'
	threads: 2
	params:
		genome="hg19",
	log:
		UNIQUEDIR+'log/{Sample}_annoCommon.log'
	shell:
		'annotatePeaks.pl {input.peaks} {params.genome} -size given\
		-annStats {output.stats} -cpu {threads} 1> {output.anno} 2> {log}'


