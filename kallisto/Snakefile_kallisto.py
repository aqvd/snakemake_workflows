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
RESDIR = config["resdir"]

if not RESDIR.endswith("kallisto/"):
	RESDIR+="kallisto/"

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
    "mm10":"" ,
    "mm38":"" ,
    "hg19":"",
    "hg38":"/data/aquevedo/resources/genomes/human/GRCh38-hg38/indexes/ensembl_103/Homo_sapiens.GRCh38.cdna.all.release-103_k31.idx",
    "-":""
    }

gtf_path = {
	"mm9" : "",
	"mm10" : "",
	"mm38" : "",
	"hg19" : "",
	"hg38" : "/data/aquevedo/resources/genomes/human/GRCh38-hg38/annotation/Homo_sapiens.GRCh38.103.gtf.gz"
	}

chr_lens : {
	"mm9" : "",
	"mm10" : "",
	"mm38" : "",
	"hg19" : "",
	"hg38" : "/data/aquevedo/resources/genomes/human/GRCh38-hg38/sequences/ensembl_103/chrLens.txt"
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

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!! Commented for inHouse data p53KO, we know the names of R1 and R2 !!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
# data["R1Path"] = [FASTQDIR + fq.replace(".fastq.gz","_R1.fastq.gz") for fq in data.File]
# data["R2Path"] = [FASTQDIR + fq.replace(".fastq.gz","_R2.fastq.gz") for fq in data.File]

data["FilePath"] = [FASTQDIR + fq for fq in data.File]

data["KallistoIx_path"] = [genome_path[i] for i in data.Genome]
data["Genome_size"] = [genome_size[i] for i in data.Genome] 
data["GTF_path"] = [gtf_path[i] for i in data.Genome] 
data["ChrLens"] = [chr_lens[i] for i in data.Genome] 

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

def get_kallisto_reads(wildcards):
	ix = data.Samples == wildcards.sample
	is_paired = data.LibraryLayout[ix].values[0] == "PAIRED"
	# If paired, return column R1 and R2 (created in line 70 and 71)
	if is_paired:
		R1_path = data.R1Path[ix].values[0]
		R2_path = data.R2Path[ix].values[0]
		return [R1_path, R2_path]
	# No paired, just return File
	else:
		# kallisto_align throw error if less than 2 files returned {input[1]} when
		# substituting the command inn the IS_PAIRED if condition
		# Return fake empty file  
		return data.FilePath[ix]

rule kallisto_align:
	input:
		get_kallisto_reads
	output:
		RESDIR + "{sample}/abundance.tsv",
		RESDIR + "{sample}/run_info.json"
	threads:
		get_resource("kallisto", "threads")
	resources:
		mem_mb = get_resource("kallisto", "mem_mb"),
		time = get_resource("kallisto", "walltime")
	params:
		is_paired = lambda wildcards: is_paired(wildcards.sample, data),
		genomeIx = lambda wildcards: expand("{genome}",genome=data.KallistoIx_path[data.Samples==wildcards.sample].values[0]),
		chrLens = lambda wildcards: data.ChrLens[data.Samples==wildcards.sample].values[0],
		reads = lambda wildcards, input: [ f"<(zcat {x})" if is_paired(wildcards.sample, data) else x for x in input ],
		fragmrent_len = 180,
		fragmrent_sd = 20,
		outdir = lambda wildcards: RESDIR + f"{wildcards.sample}/"
	conda:
		SNAKEDIR + 'envs/samtools.yaml'
	log:
		LOGDIR + "kallisto_{sample}.log"
	shell:
		'''
		if [ "{params.is_paired}" == "UNPAIRED" ]; then
		    ( kallisto quant --single --bootstrap-samples=100 \
		    	-l {params.fragment_len} -s {params.fragment_sd} \
				--index="{params.genomeIx}" --output-dir="{params.outdir}" \
				{params.reads} ) 3>&2 2>&1 1>&3 | tee {log}
		
		elif [ "{params.is_paired}" == "PAIRED" ]; then
		    ( kallisto quant --bootstrap-samples=100\
				--index="{params.genomeIx}" --output-dir="{params.outdir}" \
				{params.reads} ) 3>&2 2>&1 1>&3 | tee {log}
		fi
		'''
