import os
import  pandas as pd
import numpy as np
import glob
import re

##################################################
##				Configutation					##
##################################################
## write .yaml configuration filename
configfile: "config/config_chromHMM.yaml"


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

[tryMkdir(p) for p in (DATADIR, RESDIR, LOGDIR, RESDIR+"/fastQC")]

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
data = pd.read_csv(TABLE_NAME,sep=",")
## Add extra cols for salecting the appropriate wildcards path to files

data["Samples"] = data.Protein+"_"+data.BioProject+"_"+data.Rep.astype(str)
# Match each sample reads with appropriate input to calculate calibration factor
data["Input"] = [ data.Samples[(data.Protein=="input") & (data.BioProject==BP)].values[0] \
                 if Prot != "input" \
                 else "" \
                 for Prot,BP in zip(data.Protein,data.BioProject)  ]

data["PATH_genome"] = [genome_path[i] for i in data.Genome] 

data["Genome_size"] = [genome_size[i] for i in data.Genome]

## Remove .fastq.gz to use basename with expand() in rule "all"
data["fqBasename"] = [f.replace(".fastq.gz","") for f in data["File"]]

## function to parse config.yaml and extract resource parameters
def get_resource(rule,resource):
    try:
        return config["rules"][rule]["res"][resource]
    except KeyError:
        return config["rules"]["default"]["res"][resource]

## Decide if reference are RefSeqGenes or summits from macs2 for compute_matrix
regionsType = get_resource("compute_matrix", "regions")

if regionsType == "summits":
	regionsFile = RESDIR + '/macs/summits_unique.bed'
else:
	regionsFile = refSeq_genes_path["hg19"]

##################################################
##					RULES 						##
##################################################

## ## Rules priority 
rule all:
    input:
        ### fastQC quality control
        expand(RESDIR+"/fastQC/{fq_base}_fastqc.html", 
            fq_base=data.fqBasename.unique()),
        ## fastQC quality control
        expand(RESDIR+"/fastQC/{fq_base}_fastqc.zip", 
            fq_base=data.fqBasename.unique()),
        ## fastQ Screen contamination control
        expand(RESDIR+"/fastQScreen/{fq_base}_screen.txt", 
            fq_base=data.fqBasename.unique()),
        ## fastQ Screen contamination control
        expand(RESDIR+"/fastQScreen/{fq_base}_screen.html", 
            fq_base=data.fqBasename.unique()),
        ## .sam bowtie2 alignments. All reads aligning to reference genome
        expand(DATADIR+"/align/{sample}.sam", 
            sample=data.Samples.unique()),
        ## Macs2 peak files
        # expand(RESDIR + '/macs/{sample}_peaks.narrowPeak',
        #         sample=data.Samples[data.Protein!="input"].unique())
        ## Indexes of no dup alignments
        expand(DATADIR + "/align/{sample}_final.bai",
            sample=data.Samples.unique())


rule fastQC:
    input:
        # .fastq raw reads
        fq=FASTQDIR + '/{fq_base}.fastq.gz'
    output:
        RESDIR + "/fastQC/{fq_base}_fastqc.html",
        RESDIR + "/fastQC/{fq_base}_fastqc.zip"
    log:
        LOGDIR + "/fastQC_{fq_base}.log"
    threads:
        get_resource("fastQC", "threads")
    shell:
        "fastqc --outdir {}/fastQC --threads {{threads}} {{input.fq}} \
        |& tee {{log}}".format(RESDIR)

rule fastQScreen:
    input:
        fq=FASTQDIR + '/{fq_base}.fastq.gz'
    output:
        RESDIR+"/fastQScreen/{fq_base}_screen.txt",
        RESDIR+"/fastQScreen/{fq_base}_screen.html"
    log:
        LOGDIR + "/fastQC_{fq_base}.log"
    threads:
        get_resource("fastQScreen", "threads")
    shell:
        "fastq_screen --threads {{threads}} --aligner bowtie2 \
        --conf /home/aquevedo/opt/FastQ-Screen-0.14.1/fastq_screen.conf \
        --outdir {}/fastQScreen {{input.fq}} |& tee {{log}}".format(RESDIR)

rule bowtie2_alignTo_refGenome:
    input:
        fq=lambda wildcards: expand(FASTQDIR + '/{fq_file}', 
            fq_file=data.File[data.Samples==wildcards.sample].values)
    output:
        sam=temp(DATADIR + "/align/{sample}.sam")
    params:
        genomeIndex= lambda wildcards: expand("{genome}",
            genome=data.PATH_genome[data.Samples==wildcards.sample].values[0])
    threads:
        get_resource("bowtie2", "threads")
    resources:
        mem_mb=get_resource("bowtie2", "mem_mb")    
    log:
        LOGDIR + "/bowtie2_{sample}.log"
    run:
        reads=",".join(input.fq)
        shell('bowtie2 -x {params.genomeIndex} -U {reads} \
            -p {threads} --time --no-unal -S {output.sam} |& tee {log}')

rule sort_sam:
    input:
        sam=DATADIR + "/align/{sample}.sam"
    output:
        sorted_bam=temp(DATADIR + "/align/{sample}_sorted.bam")
        #sorted_bam=DATADIR + "/align/{sample}_sorted.bam"
    log:
        LOGDIR + "/gatk_sortSam_{sample}.log"
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
        sorted_bam=DATADIR + "/align/{sample}_sorted.bam"
    output:
        nodup_bam=DATADIR + "/align/{sample}_final.bam",
        metrics=DATADIR + "/align/stats/picardMarkDup_{sample}.txt"
    threads: 
        1
    resources:
        mem_mb=get_resource("gatk", "mem_mb"),
        walltime=get_resource("gatk","walltime")
    log:
        LOGDIR + "/gatk_markDup_{sample}.log"
    shell:
        '''
        gatk MarkDuplicates --java-options "-Xmx{resources.mem_mb}M" \
        -I {input.sorted_bam} \
        -O {output.nodup_bam} \
        -M {output.metrics} |& tee {log}
        '''

rule index_bam:
    input:
        DATADIR + "/align/{sample}_final.bam"
    output:
        DATADIR + "/align/{sample}_final.bai"
    threads: 1

    resources:
        mem_mb=get_resource("gatk", "mem_mb"),
        walltime=get_resource("gatk","walltime")
    log:
        LOGDIR + "/gatk_index_{sample}.log"
    shell:
        '''
        gatk BuildBamIndex --java-options "-Xmx{resources.mem_mb}M" \
        -I {input} |& tee {log}
        '''




