import os
import pandas as pd
import numpy as np
import glob
import re

##################################################
##				Configutation					##
##################################################
configfile: "/data/aquevedo/projects/p53_ARID1A_DKO-GSE164179/config.yaml"

FASTQDIR = config["fastqdir"]
DATADIR = config["datadir"]
RESDIR = config["resdir"]
LOGDIR = config["logdir"]

SNAKEDIR = config["snakedir"]
SCRIPTDIR=config ["scriptdir"]

TABLE_NAME = config["tableName"]

## add this to configToUse
FEATURE = config ["feature"]
ID_COUNTS = config ["id_counts"]
STRANDNESS = config ["strandness"]

##Function to create directories unless they already exist
def tryMkdir(path):
	try:
		os.makedirs(path,exist_ok=False) # raise error if exists
	except FileExistsError:
		pass

[tryMkdir(p) for p in (DATADIR, RESDIR, LOGDIR, RESDIR + "fastQC")]

## Genome hisat2 index prefixes paths
genome_path = {
	"mm9":"",
    "mm10":"/data/aquevedo/resources/genomes/mouse/mm38/GRCm38.p6_standardAndMito" ,
    "mm38":"/data/aquevedo/resources/genomes/mouse/mm38/GRCm38.p6_standardAndMito" ,
    "hg19":"",
    "hg38":"/data/aquevedo/resources/genomes/human/GRCh38-hg38/indexes/GRCh38.p12-standardAndMito/GRCh38.p12_standardAndMito",
    "-":""}

gtf_path = {
	"mm9" : "",
	"mm10" : "/data/aquevedo/resources/genomes/mouse/mm38/Mus_musculus.GRCm38.102_chrUCSC.gtf",
	"mm38" : "/data/aquevedo/resources/genomes/mouse/mm38/Mus_musculus.GRCm38.102_chrUCSC.gtf",
	"hg19" : "",
	"hg38" : "/data/aquevedo/resources/genomes/human/GRCh38-hg38/annotation/Homo_sapiens.GRCh38.97.gtf_chrUCSC.gtf"
}
## Genome sizes for big wig computation
genome_size={
	"mm9":2620345972,
    "mm10":2652783500,
    "mm38":2652783500,
    "hg19":2864785220,
    "hg38":2913022398}


##################################################
## 				READ METADATA					##
##################################################
data = pd.read_csv(TABLE_NAME,sep=",")
data[["Protein","Condition","Rep", "Run"]] = data[["Protein","Condition","Rep","Run"]].astype(str)

## Add extra cols for salecting the appropriate wildcards path to files
data["Samples"] = data.Protein +"_"+ data.Condition +"_"+ data.Rep
data["File"] = ["_".join((sample, srr)) + ".fastq.gz" for sample,srr in zip(data.Samples, data.Run)]

# raw reads filenames and R1/R2 if paired
data["fqBasename"] = [f.replace(".fastq.gz","") for f in data["File"]]
data["R1Path"] = [FASTQDIR + fq.replace(".fastq.gz","_R1.fastq.gz") for fq in data.File]
data["R2Path"] = [FASTQDIR + fq.replace(".fastq.gz","_R2.fastq.gz") for fq in data.File]
data["FilePath"] = [FASTQDIR + fq for fq in data.File]

data["HisatIx_path"] = [genome_path[i] for i in data.Genome]
data["Genome_size"] = [genome_size[i] for i in data.Genome] 
data["GTF_path"] = [gtf_path[i] for i in data.Genome] 

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
		expand(RESDIR + 'bw/{protCond}_BPM_merged.bw', 
			protCond=data.Prot_Cond[data.MergeReplicates == True].unique()),
		expand(RESDIR + 'bw/{sample}_BPM.bw', 
			sample=data.Samples[data.MergeReplicates == False].unique()),
		RESDIR + 'count_matrix.tsv'

def is_paired(sample, data):
	is_paired = data.LibraryLayout[data.Samples == sample].values[0] == "PAIRED"
	if is_paired:
		return "PAIRED"
	else:
		return "UNPAIRED"

def get_mate(fq_list, mate_id):
	for fq in fq_list:
		if fq.endswith(str(mate_id) + ".fastq.gz"):
			return fq
		else:
			return 0

def get_hisat_reads(wildcards):
	ix = data.Samples == wildcards.sample
	is_paired = data.LibraryLayout[ix].values[0] == "PAIRED"
	# If paired, return column R1 and R2 (created in line 70 and 71)
	if is_paired:
		R1_path = data.R1Path[ix].values[0]
		R2_path = data.R2Path[ix].values[0]
		return [R1_path, R2_path]
	# No paired, just return File
	else:
		# hisat_align throw error if less than 2 files returned {input[1]} when
		# substituting the command inn the IS_PAIRED if condition
		# Return fake empty file  
		return data.FilePath[ix]


rule hisat2_align:
	input:
		get_hisat_reads
	output:
		bam = DATADIR + "align/{sample}_all.bam",
		stats = DATADIR + "align/stats/{sample}_hisat2stats.txt"
	threads:
		get_resource("hisat2", "threads")
	resources:
		mem_mb = get_resource("hisat2", "mem_mb"),
		time = get_resource("hisat2", "walltime")
	params:
		is_paired = lambda wildcards: is_paired(wildcards.sample, data),
		genomeIx = lambda wildcards: expand("{genome}",genome=data.HisatIx_path[data.Samples==wildcards.sample].values[0])
	conda:
		SNAKEDIR + 'envs/samtools.yaml'
	log:
		LOGDIR + "hisat2_{sample}.log"
	shell:
		'''
		function join_by {{ local d=${{1-}} f=${{2-}}; if shift 2; then printf %s "$f" "${{@/#/$d}}"; fi; }}
		
		if [ "{params.is_paired}" == "UNPAIRED" ]; then
			U_reads=$(join_by " -U " {input})
		    
		    ( hisat2 --time --summary-file {output.stats} \
		    	--threads {threads} \
		    	-x {params.genomeIx} \
		    	-U ${{U_reads}} | \
		    samtools view -h -@ 3 -F 4 -O bam - > {output.bam} ) 3>&2 2>&1 1>&3 | \
		    tee {log}
		
		elif [ "{params.is_paired}" == "PAIRED" ]; then
			P_reads=$(join_by " -2 " {input})
		   ( hisat2 --time --summary-file {output.stats} \
			        --threads {threads} \
			        -x {params.genomeIx} \
			        -1 ${{P_reads}} | \
		    samtools view -h -@ 3 -F 12 -O bam - > {output.bam} ) 3>&2 2>&1 1>&3 | \
		    tee {log}
		fi
		'''

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
	params:
		gatk = get_resource("gatk", "bin"),
	shell:
		'''
		{params.gatk} SortSam --java-options "-Xmx{resources.mem_mb}M" \
		-I {input.bam} \
		-O {output.sorted_bam} \
		--SORT_ORDER coordinate |& tee {log}
		'''

rule index_bam:
	input:
		DATADIR + "align/{sample}_sorted.bam"
	output:
		DATADIR + "align/{sample}_sorted.bai"
	threads: 1
	resources:
		mem_mb=get_resource("gatk", "mem_mb"),
		walltime=get_resource("gatk","walltime")
	log:
		LOGDIR + "gatk/index_{sample}.log"
	params:
		gatk = get_resource("gatk", "bin")
	shell:
		'''
		{params.gatk} BuildBamIndex --java-options "-Xmx{resources.mem_mb}M" \
		-I {input} |& tee {log}
		'''

rule htseq_count:
	input:
		bams = expand(DATADIR + 'align/{sample}_sorted.bam', 
			sample = data.Samples.unique()),
		baix = expand(DATADIR + 'align/{sample}_sorted.bai', 
			sample = data.Samples.unique())
	output:
		RESDIR + 'count_matrix.tsv'
	params:
		gtf = lambda wildcards: expand("{gtf}",
			gtf=gtf_path[data.Genome.values[0]]),
		type_feature = FEATURE,
		ID_in_countMatrix = ID_COUNTS,
		is_stranded = STRANDNESS,
		tmp_counts = RESDIR + 'count_matrix.tmp'
	threads: 10
	resources:
		mem_mb = 1024 * 6,
		walltime = lambda wildcards, input: 50 * len(input.bams) 
	conda:
		SNAKEDIR + 'envs/htseq.yaml'
	log:
		LOGDIR + 'htseqCount_allSamples.txt'
	shell:
		'''
		htseq-count --format=bam --order=pos --stranded={params.is_stranded} \
		--type={params.type_feature} --idattr={params.ID_in_countMatrix} \
		--mode=union --secondary-alignments=ignore --nprocesses {threads} \
		{input.bams} {params.gtf} 2> {params.tmp_counts} 3>&2 2>&1 1>&3 | tee {log} && \
        files="{input.bams}" && \
        cat <(echo $files |xargs -n 1 basename -s "_sorted.bam" | tr "\n" "\t" |sed -E -e's/^/gene\t/') \
             {params.tmp_counts} > {output} && \
        rm {params.tmp_counts}
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
		DATADIR + "align/{Prot}_{Cond}_final_merged.bai"
	params:
		I = lambda wildcards, input: Input_merge_bam(input.bam),
		gatk = get_resource("gatk", "bin")
	threads: 2
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime")
	log:
		LOGDIR + "gatk/mergeSam_{Prot}_{Cond}.log"
	# shell:
	# 	'''
	# 	{params.gatk} MergeSamFiles --java-options "-Xmx{resources.mem_mb}M" \
	# 	--SORT_ORDER coordinate --USE_THREADING true --CREATE_INDEX true \
	# 	{params.I} --OUTPUT {output[0]} |& tee {log}
	#   '''
	run:
		if len(input.bam) > 1: # more than 1 replicate
			shell('{params.gatk} MergeSamFiles --java-options "-Xmx{resources.mem_mb}M" \
				--SORT_ORDER coordinate --USE_THREADING true --CREATE_INDEX true \
	 			{params.I} --OUTPUT {output[0]} |& tee {log}')
		else: # mjust 1 replicate
			shell("touch {output}_just1replicate.info")
			shell('echo ">> Creating BAM index for {input}"')
			shell('{params.gatk} BuildBamIndex --java-options "-Xmx{resources.mem_mb}M" \
		-I {input} |& tee {log}')
			shell('echo ">>> RENAMING {input.bam} to {output[0]}" |& tee {log}')
			shell('mv {input.bam} {output[0]} |& tee -a {log}')

rule create_bigWig:
	input:
		bam=DATADIR + "align/{sample}_sorted.bam",
		bam_index=DATADIR + "align/{sample}_sorted.bai"
	output:
		bw=RESDIR + "bw/{sample}_BPM.bw"
	params:
		genomeSize= lambda wildcards: expand("{genome_size}", 
			genome_size=data.Genome_size[data.Samples==wildcards.sample].values[0])
	threads:
		get_resource("create_bigWig", "threads")
	conda:
		SNAKEDIR + "envs/deeptools.yaml"
	resources:
		mem_mb = get_resource("create_bigWig","mem_mb"),
		walltime = get_resource("create_bigwig","walltime")
	log:
		LOGDIR + "deeptols/bamCoverage_{sample}.log"
	shell:
		'''
		bamCoverage --effectiveGenomeSize {params.genomeSize} \
		--normalizeUsing BPM -p {threads} \
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
		bw=RESDIR + "bw/{prot}_{cond}_BPM_merged.bw"
	params:
		genomeSize= lambda wildcards: expand("{genome_size}", 
			genome_size=data.Genome_size[(data.Protein==wildcards.prot) &
								   (data.Condition==wildcards.cond)].values[0])
	threads:
		get_resource("create_bigWig", "threads")
	conda:
		SNAKEDIR + "envs/deeptools.yaml"
	resources:
		mem_mb = get_resource("create_bigWig","mem_mb"),
		walltime = get_resource("create_bigwig","walltime")
	log:
		LOGDIR + "deeptols/bamCoverage_{prot}_{cond}_merged.log"
	shell:
		'''
		bamCoverage --effectiveGenomeSize {params.genomeSize} \
		--normalizeUsing BPM -p {threads} \
		-b {input.bam} -o {output.bw} |& tee {log} 
		'''
