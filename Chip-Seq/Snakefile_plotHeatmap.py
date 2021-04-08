import os
import  pandas as pd
import numpy as np
import glob
import re

##################################################
##              Configutation                   ##
##################################################
## write .yaml configuration filename
configfile: "config/test_config.yaml"


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
##              READ METADATA                   ##
##################################################
data = pd.read_csv(TABLE_NAME,sep="\t")
## Add extra cols for salecting the appropriate wildcards path to files

data["Samples"] = data.Protein +"_"+data.Condition+"_"+ data.Rep 
# Match each sample reads with appropriate input to calculate calibration factor
data["Input"] = [ data.Samples[(data.Protein=="input") & (data.Condition==Cond)].values[0] \
                 if Prot != "input" \
                 else "" \
                 for Prot,Cond in zip(data.Protein,data.Condition)  ]

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
regionsType = get_resource("compute_matrix", "regions")

if regionsType == "summits":
    regionsFile = RESDIR + '/macs/summits_unique.bed'
else:
    regionsFile = refSeq_genes_path["hg19"]

##################################################
##                  RULES                       ##
##################################################
matrixLabels = data.Samples[data.Protein != 'input'].unique()

## We can select by which samples we want to cluster the heatmap
byProt = get_resource("plot_heatmap", "cluster_by_protein")
# find indexes where Protein name appears in the filename
clustBySamp= np.where(np.array([byProt in Samp for Samp in matrixLabels ]))
# where returns tuple. 1st element is array with indexes.  
clustBySamp = " ".join(clustBySamp[0].astype(str))

rule plot_heatmap:
  input:
      matrix = RESDIR + "/deeptools/matrix_{whichSamples}.gz"
  output:
      heatmap = RESDIR + "/deeptools/heatmap_{whichSamples}.gz"
  threads:
      get_resource("plot_heatmap", "threads")
  params:
      kmeans=get_resource("plot_heatmap","kmeans"),
      clustSamp=clustBySamp
  resources:
      mem_mb=get_resource("plot_heatmap", "mem_mb"),
      walltime=get_resource("plot_heatmap", "walltime")
  conda:
      "envs/deeptools.yaml"
  shell:
      'plotHeatmap --matrixFile {input.matrix} \
      --outFileName {output.heatmap} \
      --kmeans {params.kmeans} \
      --clusterUsingSamples {params.clustSamp} \
      --sortRegions descend \
      --sortUsingSamples \
      --sortUsing mean \
      --colorMap Reds \
    '