import os
import  pandas as pd
import numpy as np
import glob
import re

##################################################
##				Configutation					##
##################################################
## write .yaml configuration filename
configfile: "config/config_meanBW.yaml"


BWDIR=config["bwdir"]
LOGDIR=config["logdir"]
TABLE_NAME=config["tableName"]

## Function to create directories unless they already exist
def tryMkdir(path):
	try:
		os.makedirs(path,exist_ok=False) # raise error if exists
	except FileExistsError:
		pass

[tryMkdir(p) for p in (BWDIR, LOGDIR)]

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

data["Samples"] = data.Protein +"_"+data.Condition+"_"+ data.Rep 
# Match each sample reads with appropriate input to calculate calibration factor
data["Input"] = [ data.Samples[(data.Protein=="input") & (data.Condition==Cond)].values[0] \
                 if Prot != "input" \
                 else "" \
                 for Prot,Cond in zip(data.Protein,data.Condition)  ]

# Input for calling peaks after merging replicates
data["InputMerged"] = [ re.sub("_S[0-9]+$","", ip) if ip != ""
                       else "" for ip in data.Input]

# All different Prot_Cond prosibilities to merge replicates
data["Prot_Cond"] = ["_".join((Prot,Cond)) for Prot,Cond in zip(data.Protein,data.Condition)]

data["PATH_genome"] = [genome_path[i] for i in data.Genome] 

data["Genome_size"] = [genome_size[i] for i in data.Genome]

data["PATH_genome_cal"] = [genome_path[i] for i in data.Norm]

data["PATH_refSeq_genes"] = [refSeq_genes_path[i] for i in data.Genome]
## Remove .fastq.gz to use basename with expand() in rule "all"
data["fqBasename"] = [f.replace(".fastq.gz","") for f in data["File"]]

data["MergeBW"]=["yes" if len(data.Rep[(data.Condition==c) & (data.Protein==p)].unique()) > 1
 else "no"
 for c,p in zip(data.Condition, data.Protein)]

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
		regionsFile = BWDIR + 'unique_summits_merged.bed'
	elif mergedRegions == "no":
		regionsFile = BWDIR + 'unique_summits_NotMerged.bed'
else:
	regionsFile = refSeq_genes_path["hg19"]


##################################################
##					RULES 						##
##################################################
rule merge_bw_only:
    input:
        expand(BWDIR + "bw/{Prot_Cond}_mean.bw" ,
            Prot_Cond=data.Prot_Cond[
            (data.Protein!="input")&
            (data.MergeBW =="yes")].unique())

rule mean_bw_scaled:
    input:
    ## bw replicates have the same Protein and Condition
        bwTreat=lambda wildcards: expand(BWDIR + "bw/{Samp}_RPKM_scaled.bw",
         Samp=data.Samples[
                         (data.MergeBW == 'yes') &
                         (data.Protein == wildcards.Prot) & 
                         (data.Condition != "siC")].unique()),
        bwCont=lambda wildcards: expand(BWDIR + "bw/{Samp}_RPKM_scaled.bw",
         Samp=data.Samples[
                         (data.MergeBW == 'yes') &
                         (data.Protein == wildcards.Prot) & 
                         (data.Condition == "siC")].unique())
    output:
        BWDIR + "bw/{Prot}_{Cond}_mean.bw"
    conda:
        'envs/deeptools.yaml'
    threads: 5
    log:
        LOGDIR + "bw/{Prot}_{Cond}_mean.log"
    shell:
        '''
        bigwigCompare --opertation mean \
        --binSize 50 \
        -b1 {input.bwTreat} \
        -b2 {input.bwCont} \
        -p {threads} \
        --outTormat bigwig \
        --outFileName {output} 
        '''

