import os
import  pandas as pd
import numpy as np
import glob
import re

##################################################
##				Configutation					##
##################################################
## write .yaml configuration filename
configfile: "config/config-MoreyPaper.yaml"


## 
FASTQDIR=config["fastqdir"]
DATADIR=config["datadir"]
RESDIR=config["resdir"]
LOGDIR=config["logdir"]
TABLE_NAME=config["tableName"]

## Function to create directories unless they already exist
def tryMkdir(path):
	try:
		os.makedirs(path,exist_ok=False) # raise error if exists
	except FileExistsError:
		pass

[tryMkdir(p) for p in (DATADIR, RESDIR, LOGDIR, RESDIR + "fastQC")]

## Genome bowtie2 index prefixes paths
genome_path = {
	"mm9":"/storage/scratch01/users/dgimenezl/genomes/mouse/mm9/mm9",
    "mm10":"/storage/scratch01/users/dgimenezl/genomes/mouse/mm10/mm10" ,
    "hg19":"/storage/scratch01/users/dgimenezl/genomes/human/hg19/hg19",
    "hg38":"/storage/scratch01/users/dgimenezl/genomes/human/hg38/hg38",
    "-":""}
refSeq_genes_path = {
	"mm9" : "",
	"mm10" : "",
	"hg19" : "/storage/scratch01/users/aquevedo/genomes/human/hg19/hg19_RefSeqCuratedGenes.bed",
	"hg38" : ""
}
## Genome sizes for big wig computation
genome_size={"mm9":2620345972,
    "mm10":2652783500,
    "hg19":2864785220,
    "hg38":2913022398}

##################################################
## 				READ METADATA					##
##################################################
data = pd.read_csv(TABLE_NAME,sep="\t")
## Add extra cols for salecting the appropriate wildcards path to files

data["Samples"] = str(data.Protein)+"_"+str(data.Condition)+"_"+ str(data.Rep) 
# Match each sample reads with appropriate input to calculate calibration factor
data["Input"] = [ data.Samples[(data.Protein=="input") & (data.Condition==Cond)].values[0] \
                 if Prot != "input" \
                 else "" \
                 for Prot,Cond in zip(data.Protein,data.Condition)  ]

# Input for calling peaks after merging replicates
data["InputMerged"] = [ re.sub("_[SR][0-9]+$","", ip) if ip != ""
                       else "" for ip in data.Input]

# All different Prot_Cond prosibilities to merge replicates
data["Prot_Cond"] = ["_".join((Prot,Cond)) for Prot,Cond in zip(data.Protein,data.Condition)]

data["PATH_genome"] = [genome_path[i] for i in data.Genome] 

data["Genome_size"] = [genome_size[i] for i in data.Genome]

data["PATH_genome_cal"] = [genome_path[i] for i in data.Norm]

data["PATH_refSeq_genes"] = [refSeq_genes_path[i] for i in data.Genome]
## Remove .fastq.gz to use basename with expand() in rule "all"
data["fqBasename"] = [f.replace(".fastq.gz","") for f in data["File"]]

## function to parse config.yaml and extract resource parameters
def get_resource(rule,resource):
    try:
        return config["rules"][rule]["res"][resource]
    except KeyError:
        return config["rules"]["default"]["res"][resource]

## Decide if reference are RefSeqGenes or summits from macs2 for compute_matrix
regionsType = config['regions']
mergedRegions = config['mergedRegions'].lower()

if regionsType == "summits":
	if mergedRegions == "yes":
		regionsFile = RESDIR + 'unique_summits_merged.bed'
	elif mergedRegions == "no":
		regionsFile = RESDIR + 'unique_summits_NotMerged.bed'
else:
	regionsFile = refSeq_genes_path["hg19"]


##################################################
##					RULES 						##
##################################################

## ## Rules priority
ruleorder: bowtie2_alignTo_calGenome > bowtie2_alignTo_refGenome 
ruleorder: unique_summits_merged > unique_summits_NotMerged
rule all:
	input:
		# ## fastQC quality control
		# expand(RESDIR + "fastQC/{fq_base}_fastqc.html", 
		# 	fq_base=data.fqBasename.unique()),
		# ## fastQC quality control
		# expand(RESDIR + "fastQC/{fq_base}_fastqc.zip", 
		# 	fq_base=data.fqBasename.unique()),
		# ## fastQ Screen contamination control
		# expand(RESDIR + "fastQScreen/{fq_base}_screen.txt", 
		# 	fq_base=data.fqBasename.unique()),
		# ## fastQ Screen contamination control
		# expand(RESDIR + "fastQScreen/{fq_base}_screen.html", 
		# 	fq_base=data.fqBasename.unique()),
		# ## .sam bowtie2 alignments. All reads aligning to reference genome
		# expand(DATADIR + "align/{sample}_all.sam", 
		# 	sample=data.Samples.unique()),
		# ## .sam bowtie2 alignments. Reads unique for reference genome
		# expand(DATADIR + "align/{sample}_onlyRef.sam", 
		# 	sample=data.Samples.unique()),
		# ## bowtie2 mapping stats 
		# expand(DATADIR + "align/stats/{sample}.txt",
		# 	sample=data.Samples.unique()),
		# ## bowtie2 maping againtst calibraton genome stats
		# expand(DATADIR + "align/stats/{sample}_calibration.txt",
		# 	sample=data.Samples[ data.PATH_genome_cal != "" ].unique()),
		# ## Calibration factors. Only for chip samples with input to compare
		# expand(DATADIR + "align/{sample}_ORI_factor.txt",
		# 	sample=data.Samples[data.Input != ""].unique()),
		# ## .bam alignment files sorted and without duplicates
		# expand(DATADIR + "align/{sample}_final.bam", 
		# 	sample=data.Samples.unique()),
		# ## Indexed .bai alignments
		# expand(DATADIR + "align/{sample}_final.bai",
		# 	sample=data.Samples.unique()),
		## Merged bams 
		expand(DATADIR + "align/{Prot_Cond}_final_merged.bam",
			Prot_Cond=data.Prot_Cond.unique()),
		## .bw files RPKM normalized
		# expand(RESDIR + "bw/{sample}_RPKM.bw", 
		# 	sample=data.Samples.unique()),
		## .bw files scaled (CPM * scaleFactor)
		# expand(RESDIR + "bw/{sample}_RPKM_scaled.bw",
		# 		sample=data.Samples[(data.PATH_genome_cal!="") & 
		# 						(data.Protein!="input")].unique()),
		## Macs2 fragment prediction for merged bam
		# expand(RESDIR + 'macs/{Prot_Cond}_predictd.txt',
		# 		Prot_Cond=data.Prot_Cond[data.Protein!="input"].unique()),
		# ## Macs2 peak files merged bam
		# expand(RESDIR + 'macs/{Prot_Cond}_peaks.narrowPeak',
		# 		Prot_Cond=data.Prot_Cond[data.Protein!="input"].unique()),
		# ## Macs2 summits files from merged bams
		# expand(RESDIR + 'macs/{Prot_Cond}_summits.bed',
		# 		Prot_Cond=data.Prot_Cond[data.Protein!="input"].unique()),
		## .bed file with all unique summits among samples
		#RESDIR + 'macs/summits_unique.bed',
		## deeptools computeMatrix out files. regionsType = "summits" or "genes"
		## look right before rule all to see how regionsType is obtained.
		# RESDIR + "deeptools/matrix_all_samples_" + regionsType + ".gz",
	#RESDIR + "deeptools/regions_all_samples_TSS_sorted.bed"

rule QC_only:
	input:
		### fastQC quality control
		expand(RESDIR + "fastQC/{fq_base}_fastqc.html", 
			fq_base=data.fqBasename.unique()),
		## fastQC quality control
		expand(RESDIR + "fastQC/{fq_base}_fastqc.zip", 
			fq_base=data.fqBasename.unique()),
		## fastQ Screen contamination control
		expand(RESDIR + "fastQScreen/{fq_base}_screen.txt", 
			fq_base=data.fqBasename.unique()),
		## fastQ Screen contamination control
		expand(RESDIR + "fastQScreen/{fq_base}_screen.html", 
			fq_base=data.fqBasename.unique())

rule merge_bams_only:
	input:
		expand(DATADIR + "align/{Prot_Cond}_final_merged.bam",
			Prot_Cond=data.Prot_Cond.unique())

rule macs2_notMerged_only:
	input:
		## Macs2 fragment prediction
		expand(RESDIR + 'macs/{sample}_predictd.txt',
				sample=data.Samples[data.Protein!="input"].unique()),
		## Macs2 peak files
		expand(RESDIR + 'macs/{sample}_peaks.narrowPeak',
				sample=data.Samples[data.Protein!="input"].unique()),
		## Macs2 summits files
		expand(RESDIR + 'macs/{sample}_summits.bed',
				sample=data.Samples[data.Protein!="input"].unique()),
		## After runing script to get unique summits save runtime here
		RESDIR + 'macs/summits_NotMerged.dateRun'


rule macs2_merged_only:
	input:
		# Macs2 fragment prediction for merged bam
		expand(RESDIR + 'macs/{Prot_Cond}_merged_predictd.txt',
				Prot_Cond=data.Prot_Cond[data.Protein!="input"].unique()),
		## Macs2 peak files merged bam
		expand(RESDIR + 'macs/{Prot_Cond}_merged_peaks.narrowPeak',
				Prot_Cond=data.Prot_Cond[data.Protein!="input"].unique()),
		## Macs2 summits files from merged bams
		expand(RESDIR + 'macs/{Prot_Cond}_merged_summits.bed',
				Prot_Cond=data.Prot_Cond[data.Protein!="input"].unique()),
		## After runing script to get unique summits save runtime here
		RESDIR + 'macs/summits_merged.dateRun'

rule bw_only:
	input:
		# .bw files RPKM normalized or scaled. One per sample. Then merge if heatmap
		# shows good correlation between samples
		expand(RESDIR + "bw/{sample}_RPKM.bw", 
			sample=data.Samples.unique()),
		# .bw files scaled (CPM * scaleFactor)
		expand(RESDIR + "bw/{sample}_RPKM_scaled.bw",
				sample=data.Samples[(data.PATH_genome_cal!="") & 
								(data.Protein!="input")].unique())

rule merge_bw_only:
	input:
		expand(RESDIR + "bw/{Prot_Cond}_mean.bw" ,
			Prot_Cond=data.Prot_Cond[data.Protein!="input"].unique())

rule compute_matrix_only:
	input:
		# deeptools computeMatrix out files. regionsType = "summits" or "genes"
		# look right before rule all to see how regionsType is obtained.
		RESDIR + "deeptools/matrix_all_samples_" + regionsType + ".gz"

rule fastQC:
	input:
		# .fastq raw reads
		fq=FASTQDIR + '{fq_base}.fastq.gz'
	output:
		RESDIR + "fastQC/{fq_base}_fastqc.html",
		RESDIR + "fastQC/{fq_base}_fastqc.zip"
	log:
		LOGDIR + "fastQC_{fq_base}.log"
	threads:
		get_resource("fastQC", "threads")
	shell:
		"fastqc --outdir {}/fastQC --threads {{threads}} {{input.fq}} \
		|& tee {{log}}".format(RESDIR)

rule fastQScreen:
	input:
		fq=FASTQDIR + '{fq_base}.fastq.gz'
	output:
		RESDIR + "fastQScreen/{fq_base}_screen.txt",
		RESDIR + "fastQScreen/{fq_base}_screen.html"
	log:
		LOGDIR + "fastQC_{fq_base}.log"
	threads:
		get_resource("fastQScreen", "threads")
	shell:
		"fastq_screen --threads {{threads}} --aligner bowtie2 \
		--conf /home/aquevedo/opt/FastQ-Screen-0.14.1/fastq_screen.conf \
		--outdir {}/fastQScreen {{input.fq}} |& tee {{log}}".format(RESDIR)
		

rule bowtie2_alignTo_refGenome:
	input:
		fq=lambda wildcards: expand(FASTQDIR + '{fq_file}', 
			fq_file=data.File[data.Samples==wildcards.sample].values)
	output:
		sam=DATADIR + "align/{sample}_onlyRef.sam",
		unal=DATADIR + "align/{sample}_unal.fastq.gz",
		stats=DATADIR + "align/stats/{sample}.txt"
	params:
		## with logical indexing retrieve the same PATH_genome n times, 
		## get the first.
		genomeIndex= lambda wildcards: expand("{genome}",
			genome=data.PATH_genome[data.Samples==wildcards.sample].values[0]),
		calGenIx = lambda wildcards: expand("{calGenome}",
			calGenome=data.PATH_genome_cal[data.Samples==wildcards.sample].values[0])
	threads:
		get_resource("bowtie2", "threads")
	resources:
		mem_mb=get_resource("bowtie2", "mem_mb")	
	log:
		LOGDIR + "bowtie2_{sample}.log"
	run:
		reads=",".join(input.fq)
		## {output.unal}: reads that do NOT align to calibration genome
		## -S /dev/null to discard aligned reads
		shell("bowtie2 -x {params.calGenIx} -U {reads} \
			-p {threads} --time	--un-gz {output.unal} -S /dev/null")
		## -S {output.sam}: reads unique to reference genome.
		## {output.stats}: align stats of reads unique to human
		shell("bowtie2 -x {params.genomeIndex} -U {output.unal} \
			-p {threads} --time	--no-unal -S {output.sam} \
			|& tee {output.stats}")

rule bowtie2_alignTo_calGenome:
	input:
		fq=lambda wildcards: expand(FASTQDIR + '{fq_file}', 
			fq_file=data.File[data.Samples==wildcards.sample].values)
	output:
		sam=DATADIR + "align/{sample}_all.sam",
		unal=DATADIR + "align/{sample}_calibration_unal.sam",
		stats= DATADIR + "align/stats/{sample}_calibration.txt"
	params:
		calGenIx = lambda wildcards: expand("{calGenome}",
			calGenome=data.PATH_genome_cal[data.Samples==wildcards.sample].values[0]),
		genomeIndex= lambda wildcards: expand("{genome}",
			genome=data.PATH_genome[data.Samples==wildcards.sample].values[0])
	threads:
		get_resource("bowtie2", "threads")
	resources:
		mem_mb=get_resource("bowtie2", "mem_mb")	
	log:
		LOGDIR + "bowtie2_calibration_{sample}.log"
	run:
		reads=",".join(input.fq)
		## Get all reads that align to reference genome in {output.sam}
		## Get reads that do NOT align to rederence in {output.unal}.
		shell("bowtie2 -x {params.genomeIndex} -U {reads} \
			-p {threads} --time	--un-gz {output.unal} -S {output.sam}")
		## {output.stats}: alignemt stats reads unique to calibration genome
		shell("bowtie2 -x {params.calGenIx} -U {output.unal} \
			-p {threads} --time	--no-unal -S /dev/null \
			|& tee {output.stats}")

rule calculate_scaled:
	input:
		# stats of alignment for Chip samples
		stats=DATADIR + "align/stats/{sample}.txt",
		calStats=DATADIR + "align/stats/{sample}_calibration.txt",
		# stats of alignment for Input samples
		inputStats=lambda wildcards: expand(
			DATADIR + "align/stats/{input_stats}.txt", 
			input_stats=data.Input[data.Samples==wildcards.sample].values[0]),
		inputCalStats=lambda wildcards: expand(
			DATADIR + "align/stats/{input_stats}_calibration.txt",
			input_stats=data.Input[data.Samples==wildcards.sample].values[0])
	output:
		DATADIR + "align/{sample}_ORI_factor.txt"
	shell:
		"n1=$(awk -f scripts/get_aligned_reads.awk {input.stats});"
		"n2=$(awk -f scripts/get_aligned_reads.awk {input.calStats});"
		"n3=$(awk -f scripts/get_aligned_reads.awk {input.inputStats});"
		"n4=$(awk -f scripts/get_aligned_reads.awk {input.inputCalStats});"
		'n=$(python -c "print(($n4 * $n1) / ($n2 * $n3))");'
		"echo $n > {output}"		

rule sort_sam:
	input:
		sam=DATADIR + "align/{sample}_all.sam"
	output:
		sorted_bam=DATADIR + "align/{sample}_sorted.bam"
		#sorted_bam=DATADIR + "align/{sample}_sorted.bam"
	log:
		LOGDIR + "gatk/sortSam_{sample}.log"
	threads: 
		1
	resources:
		mem_mb=get_resource("gatk","mem_mb"),
		walltime=get_resource("gatk","walltime")
	shell:
		'''
		gatk SortSam --java-options "-Xmx{resources.mem_mb}M" \
		-I {input.sam} \
		-O {output.sorted_bam} \
		--SORT_ORDER coordinate |& tee {log}
		'''

rule remove_duplicates:
	input:
		sorted_bam=DATADIR + "align/{sample}_sorted.bam"
	output:
		nodup_bam=DATADIR + "align/{sample}_final.bam",
		metrics=DATADIR + "align/stats/picardMarkDup_{sample}.txt"
	threads: 
		1
	resources:
		mem_mb=get_resource("gatk", "mem_mb"),
		walltime=get_resource("gatk","walltime")
	log:
		LOGDIR + "gatk/markDup_{sample}.log"
	shell:
		'''
		gatk MarkDuplicates --java-options "-Xmx{resources.mem_mb}M" \
		-I {input.sorted_bam} \
		-O {output.nodup_bam} \
		-M {output.metrics} |& tee {log}
		'''

rule index_bam:
	input:
		DATADIR + "align/{sample}_final.bam"
	output:
		DATADIR + "align/{sample}_final.bai"
	threads: 1
	resources:
		mem_mb=get_resource("gatk", "mem_mb"),
		walltime=get_resource("gatk","walltime")
	log:
		LOGDIR + "gatk/index_{sample}.log"
	shell:
		'''
		gatk BuildBamIndex --java-options "-Xmx{resources.mem_mb}M" \
		-I {input} |& tee {log}
		'''

def Input_merge_bam(bams):
    '''
    Generate formated string for -I option for MergeSamFiles
    '''
    input_bams=["-I " + str(b) for b in bams]
    input_bams=" ".join(input_bams)

    return input_bams

rule merge_bam:
	input:
		## sample replicates have same Protein and Condition
		bam=lambda wildcards: expand(DATADIR + "align/{Samp}_final.bam",
		 Samp=data.Samples[(data.Protein == wildcards.Prot) & 
		 (data.Condition == wildcards.Cond)].unique())
	output:
		## Get Prot and Cond wildcards from expanded filename  
		DATADIR + "align/{Prot}_{Cond}_final_merged.bam",
	params:
		I=lambda wildcards, input: Input_merge_bam(input.bam)
	threads: 2
	resources:
		mem_mb=get_resource("gatk", "mem_mb"),
		walltime=get_resource("gatk","walltime")
	log:
		LOGDIR + "gatk/mergeSam_{Prot}_{Cond}.log"
	shell:
		'''
		gatk MergeSamFiles --java-options "-Xmx{resources.mem_mb}M" \
		--SORT_ORDER coordinate --USE_THREADING true --CREATE_INDEX true \
		{params.I} --OUTPUT {output} |& tee {log}
		'''


rule create_bigWig:
	input:
		nodup_bam=DATADIR + "align/{sample}_final.bam",
		bam_index=DATADIR + "align/{sample}_final.bai"
	output:
		bw=RESDIR + "bw/{sample}_RPKM.bw"
	params:
		genomeSize= lambda wildcards: expand("{genome_size}", 
			genome_size=data.Genome_size[data.Samples==wildcards.sample].values[0])
	threads:
		get_resource("create_bigWig", "threads")
	conda:
		"envs/deeptools.yaml"
	resources:
		mem_mb=get_resource("create_bigWig","mem_mb"),
		walltime=get_resource("create_bigwig","walltime")
	log:
		LOGDIR + "deeptols/bamCoverage_{sample}.log"
	shell:
		'''
		bamCoverage --effectiveGenomeSize {params.genomeSize} \
		--normalizeUsing RPKM -p {threads} \
		-b {input.nodup_bam} -o {output.bw} |& tee -a {log} 
		'''

rule create_bigWig_scaled:
	input:
		nodup_bam=DATADIR + "align/{sample}_final.bam",
		bam_index=DATADIR + "align/{sample}_final.bai",
		scaleFactor=DATADIR + "align/{sample}_ORI_factor.txt"
	output:
		bw=RESDIR + "bw/{sample}_RPKM_scaled.bw"
	params:
		genomeSize= lambda wildcards: expand("{genome_size}", 
			genome_size=data.Genome_size[data.Samples==wildcards.sample].values[0])
	threads:
		get_resource("create_bigWig", "threads")
	conda:
		"envs/deeptools.yaml"
	resources:
		mem_mb=get_resource("create_bigWig","mem_mb"),
		walltime=get_resource("create_bigwig","walltime")
	log:
		LOGDIR + "deeptols/bamCoverage_{sample}_scaled.log"
	shell:
		'''
		n=$(cat {input.scaleFactor});
		bamCoverage --effectiveGenomeSize {params.genomeSize} \
		--normalizeUsing CPM --scaleFactor $n -p {threads} \
		-b {input.nodup_bam} -o {output.bw} |& tee -a {log} 
		'''		
rule merge_bw_scaled:
	input:
	## bw replicates have the same Protein and Condition
		bw=lambda wildcards: expand(RESDIR + "bw/{Samp}_RPKM_scaled.bw",
		 Samp=data.Samples[(data.Protein == wildcards.Prot) & 
		 (data.Condition == wildcards.Cond)].unique())
	output:
		RESDIR + "bw/{Prot}_{Cond}_mean.bw"
	params:
		bws = lambda wildcards, input: " ".join(input.bw)
	conda:
		'envs/wiggletools.yaml'
	threads: 1
	log:
		LOGDIR + "bw/{Prot}_{Cond}_mean.log"
	shell:
		'''
		scripts/bwMean.sh {params.bws} {output} |& tee {log}
		'''


rule macs2_notMerged_callpeak:
	input:
		treatBam= DATADIR + 'align/{sample}_final.bam',
		inputBam= lambda wildcards: expand(
					DATADIR + 'align/{inputSample}_final.bam',
			inputSample=data.Input[data.Samples == wildcards.sample].values[0])
	output:
		pred=RESDIR + 'macs/{sample}_predictd.txt',
		peaks=RESDIR + 'macs/{sample}_peaks.narrowPeak',
		summits=RESDIR + 'macs/{sample}_summits.bed'
	threads:
		get_resource("macs2", "threads")
	resources:
		mem_mb=get_resource("macs2", "mem_mb"),
		walltime=get_resource("macs2", "walltime")
	params:
		outDir=RESDIR + 'macs',
		species=get_resource("macs2", "species")
	log:
		LOGDIR + "macs/{sample}.log"
	shell:
		'scripts/macs2_callPeaks.sh {input.treatBam} {input.inputBam} \
		{params.species} {wildcards.sample} {params.outDir} &> {log}'

rule macs2_merged_callpeak:
	input:
		treatBam= lambda wildcards: expand(
			DATADIR + 'align/{Prot_Cond}_final_merged.bam',
				Prot_Cond=wildcards.Prot_Cond),
		inputBam= lambda wildcards: expand(
			DATADIR + 'align/{inputSample}_final_merged.bam',
				inputSample=data.InputMerged[data.Prot_Cond == wildcards.Prot_Cond].values[0])
	output:
		pred=RESDIR + 'macs/{Prot_Cond}_merged_predictd.txt',
		peaks=RESDIR + 'macs/{Prot_Cond}_merged_peaks.narrowPeak',
		summits=RESDIR + 'macs/{Prot_Cond}_merged_summits.bed'
	threads:
		get_resource("macs2", "threads")
	resources:
		mem_mb=get_resource("macs2", "mem_mb"),
		walltime=get_resource("macs2", "walltime")
	params:
		fileName=lambda wildcards: str(wildcards.Prot_Cond) + "_merged",
		outDir=RESDIR + 'macs',
		species=get_resource("macs2", "species")
	log:
		LOGDIR + "macs/{Prot_Cond}.log"
	shell:
		'scripts/macs2_callPeaks.sh {input.treatBam} {input.inputBam} \
		{params.species} {params.fileName} {params.outDir} &> {log}'

rule unique_summits_NotMerged:
	input:
		summits=expand(RESDIR + 'macs/{sample}_summits.bed',
			   sample=data.Samples[data.Protein != 'input'].unique())
	output:
		RESDIR + 'unique_summits_NotMerged.bed'
	threads: 1
	resources:
		mem_mb=500,
		walltime=5
	params:
		summitsBed=lambda wildcards,input: " ".join(input.summits),
		outFilename = "unique_summits_NotMerged.bed"
	log:
		LOGDIR + "macs/unique_summits.log"
	shell:
		'scripts/all_unique_summits.sh {params.outFilename} \
		{params.summitsBed} |& tee {log}'

rule unique_summits_merged:
	input:
		summits=expand(RESDIR + 'macs/{Prot_Cond}_merged_summits.bed',
			   Prot_Cond=data.Prot_Cond[data.Protein != 'input'].unique())
	output:
		RESDIR + 'unique_summits_merged.bed'
	threads: 1
	resources:
		mem_mb=500,
		walltime=5
	params:
		summitsBed=lambda wildcards,input: " ".join(input.summits),
		outFilename = "unique_summits_merged.bed"
	log:
		LOGDIR + "macs/unique_summits.log"
	shell:
		'scripts/all_unique_summits.sh {params.outFilename} \
		{params.summitsBed} |& tee {log}'



## For rule compute matrix we need to generate lists with filenames and labels
matrixFiles = data.Samples[data.Protein != 'input'].unique()
matrixLabels = " ".join(matrixFiles)

rule compute_matrix:
	input:
		bw = expand(RESDIR + "bw/{MatFile}_RPKM_scaled.bw",
				  MatFile=matrixFiles),
		regions = regionsFile ## View right before rule all
	output:
		matrix = RESDIR + "deeptools/matrix_all_samples_"+regionsType+".gz",
		#sortRegions = RESDIR + "deeptools/regions_all_samples_TSS_sorted.bed"
	threads:
		get_resource("compute_matrix", "threads")
	params:
		labels = matrixLabels
	conda:
		"envs/deeptools.yaml"
	resources:
		mem_mb=get_resource("compute_matrix", "mem_mb"),
		walltime=get_resource("compute_matrix", "walltime")
	shell:
		'computeMatrix reference-point --referencePoint center \
			--scoreFileName {input.bw} \
			--regionsFileName {input.regions} \
			--upstream 2500 --downstream 2500 \
			--binSize 50 \
			--missingDataAsZero \
			--samplesLabel {params.labels} \
			--numberOfProcessors {threads} \
			--outFileName {output.matrix} '


