import os
import  pandas as pd
import numpy as np
import glob
import re

##################################################
##				Configutation					##
##################################################
configfile: "config/config-ToUse.yaml"

FASTQDIR = config["fastqdir"]
DATADIR = config["datadir"]
RESDIR = config["resdir"]
LOGDIR = config["logdir"]
TABLE_NAME = config["tableName"]
## add this to configToUse
GTF = config["gtf"]
FEATURE = config ["feature"]
ID_COUNTS = config ["id_counts"]
STRANDNESS = config ["strandness"]

## Function to create directories unless they already exist
def tryMkdir(path):
	try:
		os.makedirs(path,exist_ok=False) # raise error if exists
	except FileExistsError:
		pass

[tryMkdir(p) for p in (DATADIR, RESDIR, LOGDIR, RESDIR + "fastQC")]

## Genome bowtie2 index prefixes paths
genome_path = {
	"mm9":"",
    "mm10":"" ,
    "hg19":"/storage/scratch01/users/aquevedo/genomes/human/hg19/hisat/hg19/genome",
    "hg38":"",
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

data["Samples"] = data.Protein +"_"+ data.Condition +"_"+ data.Rep

data["fqBasename"] = [f.replace(".fastq.gz","") for f in data["File"]]

data["HisatIx_path"] = [genome_path[i] for i in data.Genome]
data["Genome_size"] = [genome_size[i] for i in data.Genome] 

# All different Prot_Cond prosibilities to merge replicates
data["Prot_Cond"] = ["_".join((Prot,Cond)) for Prot,Cond in zip(data.Protein,data.Condition)]

# Column of boolean to decide if merge replicates
def multipleReps(sample):
    prot, cond, _ = sample.split("_")
    reps = len(data.Rep[(data.Protein == prot) & (data.Condition == cond)].unique())
    
    return reps > 1

data["MergeReplicates"]=[multipleReps(sample) for sample in data.Samples]

## Regions Used for ComputeMatrix
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


def get_resource(rule,resource):
    try:
        return config["rules"][rule]["res"][resource]
    except KeyError:
        return config["rules"]["default"]["res"][resource]

##################################################
##					RULES 						##
##################################################
rule all:
	input:
		expand(RESDIR + 'bw/{protCond}_RPKM_merged.bw', 
			protCond=data.Prot_Cond[data.MergeReplicates == True].unique()),
		expand(RESDIR + 'bw/{sample}_RPKM.bw', 
			sample=data.Samples[data.MergeReplicates == False].unique()),
		RESDIR + 'count_matrix.tsv'

rule hisat2_align_and_sortBam:
	input:
		fq = lambda wildcards: expand(FASTQDIR + '{fq_file}', 
			fq_file=data.File[data.Samples==wildcards.sample].values)
	output:
		bam = DATADIR + "align/{sample}_sorted.bam",
		bai = DATADIR + "align/{sample}_sorted.bam.bai",
		stats = DATADIR + "align/stats/{sample}_hisat2stats.txt"
	params:
		genomeIx = lambda wildcards: expand("{genome}",
			genome=data.HisatIx_path[data.Samples==wildcards.sample].values[0])
	threads: 
		get_resource("hisat2", "threads")
	resources:
		mem_mb = get_resource("hisat2", "mem_mb")
	log:
		LOGDIR + "hisat2_{sample}.log"
	run:
		reads = ','.join(input.fq)
		shell('\
			(hisat2 --time --summary-file {output.stats} \
			--no-unal --threads {threads} \
			-x {params.genomeIx} \
			-U ' + reads + ' \
			-S - | \
			samtools sort -@ 6 -O bam - > {output.bam} ) 3>&2 2>&1 1>&3 | \
			tee {log}')
		shell('samtools index {output.bam} |& tee -a {log}')

rule htseq_count:
	input:
		bams = expand(DATADIR + 'align/{sample}_sorted.bam', 
			sample = data.Samples.unique())
	output:
		RESDIR + 'count_matrix.tsv'
	params:
		gtf = GTF ,
		type_feature = FEATURE,
		ID_in_countMatrix = ID_COUNTS,
		is_stranded = STRANDNESS
	threads: 1
	resources:
		mem_mb = 2000,
		walltime = lambda wildcards, input: 40 * len(input.bams) 
	conda:
		'envs/htseq.yaml'
	log:
		LOGDIR + 'htseqCount_allSamples.txt'
	shell:
		'''
		htseq-count --format=bam --order=pos --stranded={params.is_stranded} \
		--type={params.type_feature} --idattr={params.ID_in_countMatrix} \
		--mode=union --secondary-alignments=ignore \
		{input.bams} {params.gtf} > {output} 3>&1 2>&1 1>&3 | tee {log}
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
		bam=lambda wildcards: expand(DATADIR + "align/{Samp}_sorted.bam",
			Samp=data.Samples[(data.Protein == wildcards.Prot) & 
				(data.Condition == wildcards.Cond)].unique())
	output:
		## Get Prot and Cond wildcards from expanded filename  
		DATADIR + "align/{Prot}_{Cond}_final_merged.bam",
	params:
		I = lambda wildcards, input: Input_merge_bam(input.bam)
	threads: 2
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime")
	log:
		LOGDIR + "gatk/mergeSam_{Prot}_{Cond}.log"
	# shell:
	# 	'''
	# 	gatk MergeSamFiles --java-options "-Xmx{resources.mem_mb}M" \
	# 	--SORT_ORDER coordinate --USE_THREADING true --CREATE_INDEX true \
	# 	{params.I} --OUTPUT {output} |& tee {log}
	#   '''
	run:
		if len(input.bam) > 1: # more than 1 replicate
			shell('gatk MergeSamFiles --java-options "-Xmx{resources.mem_mb}M" \
				--SORT_ORDER coordinate --USE_THREADING true --CREATE_INDEX true \
	 			{params.I} --OUTPUT {output} |& tee {log}')
		else: # mjust 1 replicate
			shell("touch {output}_just1replicate.info")
			shell('echo ">> Creating BAM index for {input}"')
			shell('gatk BuildBamIndex --java-options "-Xmx{resources.mem_mb}M" \
		-I {input} |& tee {log}')
			shell('echo ">>> RENAMING {input.bam} to {output}" |& tee {log}')
			shell('mv {input.bam} {output} |& tee -a {log}')

rule create_bigWig:
	input:
		bam=DATADIR + "align/{sample}_sorted.bam",
		bam_index=DATADIR + "align/{sample}_sorted.bam.bai"
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
		mem_mb = get_resource("create_bigWig","mem_mb"),
		walltime = get_resource("create_bigwig","walltime")
	log:
		LOGDIR + "deeptols/bamCoverage_{sample}.log"
	shell:
		'''
		bamCoverage --effectiveGenomeSize {params.genomeSize} \
		--normalizeUsing RPKM -p {threads} \
		-b {input.bam} -o {output.bw} |& tee {log} 
		'''

rule create_bigWig_mergedReps:
	input:
		bam=lambda wildcards: expand(DATADIR + "align/{ProtCond}_final_merged.bam",
			ProtCond=data.Prot_Cond[(data.Protein==wildcards.prot) &
								   (data.Condition==wildcards.cond)].values[0]),
		bam_index=lambda wildcards: expand(DATADIR + "align/{ProtCond}_final_merged.bai",
			ProtCond=data.Prot_Cond[(data.Protein==wildcards.prot) &
								   (data.Condition==wildcards.cond)].values[0]),
	output:
		bw=RESDIR + "bw/{prot}_{cond}_RPKM_merged.bw"
	params:
		genomeSize= lambda wildcards: expand("{genome_size}", 
			genome_size=data.Genome_size[data.Samples==wildcards.sample].values[0])
	threads:
		get_resource("create_bigWig", "threads")
	conda:
		"envs/deeptools.yaml"
	resources:
		mem_mb = get_resource("create_bigWig","mem_mb"),
		walltime = get_resource("create_bigwig","walltime")
	log:
		LOGDIR + "deeptols/bamCoverage_{prot}_{cond}_merged.log"
	shell:
		'''
		bamCoverage --effectiveGenomeSize {params.genomeSize} \
		--normalizeUsing RPKM -p {threads} \
		-b {input.bam} -o {output.bw} |& tee {log} 
		'''
