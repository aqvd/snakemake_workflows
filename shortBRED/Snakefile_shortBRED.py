import os
import pandas as pd
import numpy as np
import glob
import re

##################################################
##				Configutation					##
##################################################
## write .yaml configuration filename
configfile: "/work_beegfs/sukmb552/snakemake_workflows/shortBRED/config.yaml"

TABLE_NAME = config["tableName"]

# Outpur directories
FASTQDIR = config["fastqdir"]
DATADIR = config["datadir"]
RESDIR = config["resdir"]
LOGDIR = config["logdir"]

# Software directories
BBTOOLS_FOLDER = config["bbtools_folder"]

# Dir with Conda ennvironment .yaml definintions and tmporary folder
TMP_FOLDER = config["tmp_folder"]
#SCRIPT_FOLDER = config["script_folder"]

## Function to create directories unless they already exist
def tryMkdir(path):
	try:
		os.makedirs(path,exist_ok=False) # raise error if exists
	except FileExistsError:
		pass

[tryMkdir(p) for p in (DATADIR + "stats", LOGDIR, TMP_FOLDER,
					   RESDIR + "fastQC")]

##################################################
## 				    METADATA					##
##################################################
def rep_undersc(col, rep = ":"):
	return [str(i).replace("_", rep) for i in col.tolist()]

data = pd.read_csv(TABLE_NAME,sep="\t")

## R1.fasq.gz filename maped to Sample metadata is what we require from metadata
## mate2 we will repalce from gf
data["fqR2"] = [re.sub(r"_R1_001","_R2_001", f) for f in data["fqR1"]]

data[["SampleName1","SampleName2","SampleName3"]] = data[["SampleName1","SampleName2","SampleName3"]].apply(rep_undersc, axis = 0)

data["SampleName"] = data.SampleName1+"-"+ data.SampleName2 + "-"+ data.SampleName3

# save metadata used
metadata_saved = os.path.join(DATADIR, "metadata_used_by_snakemake.csv")
data.to_csv(metadata_saved, sep = "\t", index = False)

# chenck unique sample names
_, c = np.unique(data["SampleName"], return_counts=True)

assert all(c==1), f"Duplicated Samples!\nCheck Snamekemale used metadata file: {metadata_saved}" 


##################################################
##					RESOURCES					##
##################################################

## function to parse config.yaml and extract resource parameters
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
		expand(RESDIR + "{sample}_tmp/markers.tab", sample = data.SampleName.unique()),


def get_readPair(pairID, fq_list):
	"""
	In case several .fastq per read mate, get all to concatenate with zcat prior bwa
	"""
	search = [re.match(r".+_{}_.+".format(pairID), fq) for fq in fq_list]
	filtered = list(filter(lambda x: x != None, search))

	return " ".join([filt[0] for filt in filtered])
	

rule bbduk_trim_adapt:
	input:
		R1 = lambda wildcards: expand('{R1}',
			R1 = data.fqR1[data.SampleName == wildcards.sample].values),
		R2 = lambda wildcards: expand('{R2}',
			R2 = data.fqR2[data.SampleName == wildcards.sample].values)
	output:
		trim1 = temp( DATADIR + "{sample}_trim_R1.fq.gz"),
		trim2 = temp( DATADIR + "{sample}_trim_R2.fq.gz")
	log:
		LOGDIR + '{sample}_bbduk_trim.log'
	threads:
		get_resource("bbduk_trim", "threads")
	resources:
		mem_mb = get_resource("bbduk_trim", "mem_mb"),
		walltime = get_resource("bbduk_trim", "walltime")
	params:
		adapters=BBTOOLS_FOLDER + "/resources/adapters.fa",
		bbt_folder = BBTOOLS_FOLDER
	shell:
		'''
		{params.bbt_folder}/bbduk.sh -Xmx{resources.mem_mb}m \
        in1={input.R1} in2={input.R2} \
        out1={output.trim1} out2={output.trim2} \
        ref={params.adapters} \
        ktrim=r k=23 mink=11 hdist=1 tpe tbo \
        qtrim=r trimq=10 \
        maq=10 \
        pigz=t unpigz=t |& tee {log}
		'''

rule bbduk_phiX:
	input:
		R1 =  DATADIR + "{sample}_trim_R1.fq.gz",
		R2 =  DATADIR + "{sample}_trim_R2.fq.gz"
	output:
		trim1 = temp( DATADIR + "{sample}_trim_noPhix_R1.fq.gz"),
		trim2 = temp( DATADIR + "{sample}_trim_noPhix_R2.fq.gz")
	log:
		LOGDIR + '{sample}_bbduk_phiX.log'
	threads:
		get_resource("bbduk_trim", "threads")
	resources:
		mem_mb = get_resource("bbduk_trim", "mem_mb"),
		walltime = get_resource("bbduk_trim", "walltime")
	params:
		bbt_folder = BBTOOLS_FOLDER,
		phiX= BBTOOLS_FOLDER + "/resources/phix174_ill.ref.fa.gz",
		stats= LOGDIR + '{sample}_bbduk_phiX_stats.txt'
	shell:
		'''
		{params.bbt_folder}/bbduk.sh -Xmx{resources.mem_mb}m \
        in1={input.R1} in2={input.R2} \
        out1={output.trim1} out2={output.trim2} \
        ref={params.phiX} \
        ktrim=r k=23 mink=11 hdist=1 tpe tbo \
        k=31 hdist=1 stats={params.stats} \
        pigz=t unpigz=t" |& tee {log}
		'''

rule bbmap_host:
	input:
		in1 = DATADIR + "{sample}_trim_noPhix_R1.fq.gz",
		in2 = DATADIR + "{sample}_trim_noPhix_R2.fq.gz"
	output:
		DATADIR + "{sample}_trim_noPhix_noHost_inlv.fasta.gz"
	log:
		LOGDIR + '{sample}_bbmap_host.log'
	threads:
		get_resource("bbmap", "threads")
	resources:
		mem_mb = get_resource("bbmap", "mem_mb"),
		walltime = get_resource("bbmap", "walltime")
	params:
		bbt_folder = BBTOOLS_FOLDER,
		ix = get_resource("bbmap", "human_BBmapIX_dir"),
		stats=LOGDIR + '/{sample}_bbmap_host_stats.txt'
	shell:
		'''
		{params.bbt_folder}/bbmap.sh -Xmx${resources.mem_mb}M \
        fast \
        threads={threads} \
        minratio=0.9 maxindel=3 bwr=0.16 bw=12 minhits=2 qtrim=r trimq=10 \
        untrim idtag printunmappedcount kfilter=25 maxsites=1 k=14 \
        in1={input.in1} \
        in2={input.in2} \
        outu={output} \
        outm=/dev/null \
        path={params.ix} \
        pigz=t unpigz=t" |& tee {log}
		'''


rule shortBRED_quant:
	input:
		fa = DATADIR + "{sample}_trim_noPhix_noHost_inlv.fasta.gz"
	output:
		counts = RESDIR + "{sample}_tmp/markers.tab"
	threads:
		get_resource("shortBRED_quant", "threads")
	resources:
		mem_mb = get_resource("bbduk_map", "mem_mb"),
		walltime = get_resource("bbduk_map", "walltime")
	params:
		markers=get_resource("shortBRED_quant", "markers"),
		usearch_path=get_resource("shortBRED_quant", "usearch_path"),
		tmpdir = lambda wildcards: RESDIR + f"{wildcards.sample}_tmp"
	log:
		LOGDIR + "{sample}_shortBREDquant.log"
	shell:
		'''
		conda activate biobakery && \
		shortbred_quantify.py \
        --markers {params.markers} \
        --wgs {input.fa} \
        --usearch {params.usearch_path} \
        --threads {threads} \
        --results {output.counts} --tmp {params.tmpdir} |& tee {log}
		'''
