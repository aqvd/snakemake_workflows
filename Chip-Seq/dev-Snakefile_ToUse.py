import os
import  pandas as pd
import numpy as np
import glob
import re

##################################################
##				Configutation					##
##################################################
## write .yaml configuration filename
configfile: "config/config-ToUse.yaml"


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
    "-":"NoCalibration"}

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
## 				    METADATA					##
##################################################
data = pd.read_csv(TABLE_NAME,sep="\t")

data[["Protein","Condition","Rep"]] = data[["Protein","Condition","Rep"]].astype(str)

## Add extra cols for salecting the appropriate wildcards path to files
data["Samples"] = data.Protein +"_"+ data.Condition +"_"+ data.Rep

# Match each sample reads with appropriate input to calculate calibration factor
data["Input"] = [ data.Samples[(data.Protein=="input") & (data.Condition==Cond)].values[0] \
                 if Prot != "input" and Cond in data.Condition[data.Protein=="input"].values \
                 else "" \
                 for Prot,Cond in zip(data.Protein,data.Condition)  ]

# Input for calling peaks after merging replicates. Remove "Rep" part from Input
data["InputMerged"] = [ re.sub("^(.+)_(.+)_(.+)$",r"\1_\2", ip) if ip != ""
                       else "" for ip in data.Input]

## Create column to know which libraries have more than 1 replicate and must
## be merged
def multipleReps(sample):
    prot, cond, _ = sample.split("_")
    reps= len(data.Rep[(data.Protein==prot) & (data.Condition==cond)].unique())
    
    return reps > 1

data["MergeReplicates"] = [multipleReps(sample) for sample in data.Samples]
# All different Prot_Cond prosibilities to merge replicates

data["Prot_Cond"] = ["_".join((Prot,Cond)) for Prot,Cond in zip(data.Protein,data.Condition)]

data["PATH_genome"] = [genome_path[i] for i in data.Genome] 

data["Genome_size"] = [genome_size[i] for i in data.Genome]

data["PATH_genome_cal"] = [genome_path[i] for i in data.Norm]

data["PATH_refSeq_genes"] = [refSeq_genes_path[i] for i in data.Genome]

## Remove .fastq.gz to use basename with expand() in rule "all"
data["fqBasename"] = [f.replace(".fastq.gz","") for f in data["File"]]

##################################################
##					RESOURCES					##
##################################################

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

print(data.Samples.unique())

##################################################
##					RULES 						##
##################################################

## ## Rules priority
ruleorder: bowtie2_alignTo_calGenome > bowtie2_alignTo_refGenome 
ruleorder: unique_summits_merged > unique_summits_NotMerged
rule all:
	input:
		# ## .sam bowtie2 alignments. All reads aligning to reference genome
		# expand(DATADIR + "align/{sample}_all.bam", 
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
		# fastQ Screen contamination control
		expand(RESDIR + "fastQScreen/{fq_base}_screen.txt", 
			fq_base=data.fqBasename.unique()),
		## fastQ Screen contamination control
		expand(RESDIR + "fastQScreen/{fq_base}_screen.html", 
			fq_base=data.fqBasename.unique())

rule noDup_alignments:
	input:
		expand(DATADIR + "align/{sample}_final.bam", 
		 	sample=data.Samples.unique())

rule merge_bams_only:
	input:
		expand(DATADIR + "align/{Prot_Cond}_final_merged.bam",
			Prot_Cond=data.Prot_Cond[data.MergeReplicates == True].unique())

rule macs2Callpeak_and_BigWig:
	input:
	## No need to merge replicates
		## Macs2 fragment prediction
		expand(RESDIR + 'macs/{sample}_predictd.txt',
				sample=data.Samples[(data.Protein!="input") & 
									(data.MergeReplicates==False)].unique()),
		## Macs2 peak files
		expand(RESDIR + 'macs/{sample}_peaks.narrowPeak',
				sample=data.Samples[(data.Protein!="input") & 
									(data.MergeReplicates==False)].unique()),
		## Macs2 summits files
		expand(RESDIR + 'macs/{sample}_summits.bed',
				sample=data.Samples[(data.Protein!="input") & 
									(data.MergeReplicates==False)].unique()),

	## When merge replicates is needed
		# Macs2 fragment prediction for merged bam
		expand(RESDIR + 'macs/{Prot_Cond}_merged_predictd.txt',
				Prot_Cond=data.Prot_Cond[(data.Protein!="input") & 
										(data.MergeReplicates==True)].unique()),
		## Macs2 peak files merged bam
		expand(RESDIR + 'macs/{Prot_Cond}_merged_peaks.narrowPeak',
				Prot_Cond=data.Prot_Cond[(data.Protein!="input") & 
										(data.MergeReplicates==True)].unique()),
		## Macs2 summits files from merged bams
		expand(RESDIR + 'macs/{Prot_Cond}_merged_summits.bed',
				Prot_Cond=data.Prot_Cond[(data.Protein!="input") & 
										(data.MergeReplicates==True)].unique()),

		# .bw files RPKM normalized or scaled. One per sample. Then merge if heatmap
		# shows good correlation between samples
		expand(RESDIR + "bw/{sample}_RPKM.bw", 
			sample=data.Samples.unique()),
		#.bw files scaled (CPM * scaleFactor)
		expand(RESDIR + "bw/{sample}_RPKM_scaled.bw",
				sample=data.Samples[(data.PATH_genome_cal!="NoCalibration") & 
								(data.Protein!="input")].unique()),
		## .bw files merging replicates when there is no calibration
		expand(RESDIR + "bw/{Prot_Cond}_RPKM_merged.bw",
				Prot_Cond=data.Prot_Cond[(data.PATH_genome_cal=="NoCalibration") & 
									(data.MergeReplicates==True)].unique()),

# rule merge_bw_only:
# 	input:
# 		expand(RESDIR + "bw/{Prot_Cond}_mean.bw" ,
# 			Prot_Cond=data.Prot_Cond[data.Protein!="input"].unique())

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
		"fastqc --outdir {}fastQC --threads {{threads}} {{input.fq}} \
		|& tee {{log}}".format(RESDIR)

rule fastQScreen:
	input:
		fq=FASTQDIR + '{fq_base}.fastq.gz'
	output:
		RESDIR + "fastQScreen/{fq_base}_screen.txt",
		RESDIR + "fastQScreen/{fq_base}_screen.html"
	log:
		LOGDIR + "fastQScreen_{fq_base}.log"
	threads:
		get_resource("fastQScreen", "threads")
	resources:
		mem_mb=get_resource("fastQScreen", "mem_mb")
	shell:
		"fastq_screen --threads {{threads}} --aligner bowtie2 \
		--conf /home/aquevedo/opt/FastQ-Screen-0.14.1/fastq_screen.conf \
		--outdir {}fastQScreen {{input.fq}} |& tee {{log}}".format(RESDIR)
		

rule bowtie2_alignTo_refGenome:
	input:
		fq=lambda wildcards: expand(FASTQDIR + '{fq_file}', 
			fq_file=data.File[data.Samples==wildcards.sample].values)
	output:
		stats=DATADIR + "align/stats/{sample}_onlyRef.txt"
	params:
		## temporary file
		tmp_unal=lambda wildcards: expand(DATADIR + "align/{sample}_unal.fastq.gz",
			sample=wildcards.sample),
		## with logical indexing retrieve the same PATH_genome ntimes, 
		## get the first.
		genomeIndex= lambda wildcards: expand("{genome}",
			genome=data.PATH_genome[data.Samples==wildcards.sample].values[0]),
		calGenIx = lambda wildcards: expand("{calGenome}",
			calGenome=data.PATH_genome_cal[data.Samples==wildcards.sample].values[0])
	threads:
		get_resource("bowtie2", "threads")
	resources:
		mem_mb=get_resource("bowtie2", "mem_mb"),
		walltime=get_resource("bowtie2", "walltime")
	log:
		LOGDIR + "bowtie2_{sample}.log"
	run:
		reads=",".join(input.fq)
		## {output.unal}: reads that do NOT align to calibration genome
		## -S /dev/null to discard aligned reads
		shell("bowtie2 -x {params.calGenIx} -U {reads} \
			-p {threads} --time	--un-gz {params.tmp_unal} -S /dev/null |& tee {log}")
		## -S {output.sam}: reads unique to reference genome.
		## {output.stats}: align stats of reads unique to human
		shell("bowtie2 -x {params.genomeIndex} -U {params.tmp_unal} \
			-p {threads} --time	--no-unal -S /dev/null |& \
			tee {output.stats} && \
			rm {params.tmp_unal}")

rule bowtie2_alignTo_calGenome:
	input:
		fq=lambda wildcards: expand(FASTQDIR + '{fq_file}', 
			fq_file=data.File[data.Samples==wildcards.sample].values)
	output:
		bam=DATADIR + "align/{sample}_all.bam",
		stats= DATADIR + "align/stats/{sample}_calibration.txt"
	params:
		## temporary file
		tmp_unal=lambda wildcards: expand(DATADIR + "align/{sample}_unal.fastq.gz",
			sample=wildcards.sample),
		calGenIx = lambda wildcards: expand("{calGenome}",
			calGenome=data.PATH_genome_cal[data.Samples==wildcards.sample].values[0]),
		genomeIndex= lambda wildcards: expand("{genome}",
			genome=data.PATH_genome[data.Samples==wildcards.sample].values[0]),
		reads =lambda wildcards, input: ' '.join(input.fq),
		stat_dir=DATADIR + "align/stats"
	threads:
		get_resource("bowtie2", "threads") - 3
	resources:
		mem_mb=get_resource("bowtie2", "mem_mb"),
		walltime=get_resource("bowtie2", "walltime")
	conda:
		'envs/samtools.yaml'	
	log:
		LOGDIR + "bowtie2_calibration_{sample}.log"
	shell:
		"scripts/bowtie2_alignTo_calGenome.sh {output.bam} \
		{output.stats} {params.stat_dir} {params.tmp_unal} \
		{params.calGenIx} {params.genomeIndex} \
		{threads} {log} {params.reads}"

rule calculate_scaled:
	input:
		# stats of alignment for Chip samples
		stats=DATADIR + "align/stats/{sample}_onlyRef.txt",
		calStats=DATADIR + "align/stats/{sample}_calibration.txt",
		# stats of alignment for Input samples
		inputStats=lambda wildcards: expand(
			DATADIR + "align/stats/{input_stats}_onlyRef.txt", 
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
		bam=DATADIR + "align/{sample}_all.bam"
	output:
		sorted_bam=DATADIR + "align/{sample}_sorted.bam"
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
		-I {input.bam} \
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
		--REMOVE_DUPLICATES true \
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
				(data.Condition == wildcards.Cond) & 
				(data.MergeReplicates==True)].unique())
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
	# shell:
	# 	'''
	# 	gatk MergeSamFiles --java-options "-Xmx{resources.mem_mb}M" \
	# 	--SORT_ORDER coordinate --USE_THREADING true --CREATE_INDEX true \
	# 	{params.I} --OUTPUT {output[0]} |& tee {log}
	#   '''
	run:
		if len(input.bam) > 1:
			shell('gatk MergeSamFiles --java-options "-Xmx{resources.mem_mb}M" \
				--SORT_ORDER coordinate --USE_THREADING true --CREATE_INDEX true \
	 			{params.I} --OUTPUT {output[0]} |& tee {log}')
		else:
			shell("touch {output}_just1replicate.info")
			shell('echo ">>> RENAMING {input.bam} to {output[0]}" |& tee {log}')
			shell('mv {input.bam} {output[0]} |& tee -a {log}')


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
		--ignoreDuplicates \
		--normalizeUsing CPM --scaleFactor $n -p {threads} \
		-b {input.nodup_bam} -o {output.bw} |& tee -a {log} 
		'''	

rule create_bw_mergedReplicates:
	input:
		merged_bam = lambda wildcards: expand(
				DATADIR + "align/{prot_cond}_final_merged.bam",
					prot_cond=wildcards.Prot_Cond),
	output:
		bw = RESDIR + "bw/{Prot_Cond}_RPKM_merged.bw",
	params:
		genomeSize= lambda wildcards: expand("{genome_size}", 
			genome_size=data.Genome_size[(data.Prot_Cond==wildcards.Prot_Cond)].values[0])
	threads:
		get_resource("create_bigWig", "threads")
	conda:
		"envs/deeptools.yaml"
	resources:
		mem_mb=get_resource("create_bigWig","mem_mb"),
		walltime=get_resource("create_bigwig","walltime")
	log:
		LOGDIR + "deeptols/bamCoverage_{Prot_Cond}_merged.log"
	shell:
		'''
		bamCoverage --effectiveGenomeSize {params.genomeSize} \
		--normalizeUsing RPKM -p {threads} \
		-b {input.merged_bam} -o {output.bw} |& tee -a {log} 
		'''

rule create_bigWig_InputNorm:
	input:
		nodup_bam=DATADIR + "align/{sample}_final.bam",
		bam_index=DATADIR + "align/{sample}_final.bai",
		nodup_bam_input=lambda wildcards: expand(DATADIR + "align/{input}_final.bam",
			input=data.Input[data.Samples==wildcards.sample].values[0]),
		bam_index_input=lambda wildcards: expand(DATADIR + "align/{input}_final.bai",
			input=data.Input[data.Samples==wildcards.sample].values[0])
	output:
		bw=RESDIR + "bw/{sample}_BPM_inputNormalized.bw"
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
		LOGDIR + "deeptols/bwInputNormalized_{sample}.log"
	shell:
	# BPM normalization as RPKM not good for between samples comparison
		'''
		bamCompare -b1 {input.nodup_bam} -b2 {input.nodup_bam_input} \
		-o {output.bw} --outFileFormat bigwig \
		--binSize 50 \
		--ignoreDuplicates \
		--scaleFactorsMethod None \
		--effectiveGenomeSize {params.genomeSize} \
		--normalizeUsing BPM -p {threads} |& tee -a {log} 
		'''


rule macs2_notMerged_callpeak:
	input:
		treatBam = DATADIR + 'align/{sample}_final.bam',
		treatStats = DATADIR + "align/stats/picardMarkDup_{sample}.txt",
		inputBam = lambda wildcards: expand(
			DATADIR + 'align/{inputSample}_final.bam',
			inputSample=data.Input[data.Samples == wildcards.sample].values[0]),
		inputStats = lambda wildcards: expand(
			DATADIR + "align/stats/picardMarkDup_{inputSample}.txt",
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
		species=get_resource("macs2", "species"),
		inputSample=lambda wildcards: data.Input[data.Samples == wildcards.sample].values[0]
	log:
		LOGDIR + "macs/{sample}.log"

	conda:
		"envs/bc-1_07.yaml"
	shell:
		'scripts/macs2_callPeaks_downsampling.sh {input.treatBam} {input.inputBam} \
		{params.species} {wildcards.sample} {params.inputSample} \
		{params.outDir} ' + DATADIR + ' |& tee {log}'

rule macs2_merged_callpeak:
	input:
		treatBam= lambda wildcards: expand(
			DATADIR + 'align/{Prot_Cond}_final_merged.bam',
				Prot_Cond=wildcards.Prot_Cond),
		inputBam= lambda wildcards: expand(
			DATADIR + 'align/{inputSample}_final_merged.bam',
				inputSample=data.InputMerged[(data.Prot_Cond==wildcards.Prot_Cond)&
				(data.MergeReplicates==True) ].values[0])
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
		'scripts/macs2_callPeaks_backup.sh {input.treatBam} {input.inputBam} \
		{params.species} {params.fileName} {params.outDir} |& tee {log}'

rule unique_summits_NotMerged:
	input:
		summits=expand(RESDIR + 'macs/{sample}_summits.bed',
			   sample=data.Samples[data.Protein != 'input'].unique())
	output:
		RESDIR + 'macs/unique_summits_NotMerged.bed'
	threads: 1
	resources:
		mem_mb=2000,
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
			   Prot_Cond=data.Prot_Cond[ (data.Protein != 'input') & 
				(data.MergeReplicates==True) ].unique())
	output:
		RESDIR + 'macs/unique_summits_merged.bed'
	threads: 1
	resources:
		mem_mb=2000,
		walltime=5
	params:
		summitsBed=lambda wildcards,input: " ".join(input.summits),
		outFilename = "unique_summits_merged.bed"
	log:
		LOGDIR + "macs/unique_summits.log"
	shell:
		'scripts/all_unique_summits.sh {params.outFilename} \
		{params.summitsBed} |& tee {log}'
