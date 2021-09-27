import os
import pandas as pd
import numpy as np
import glob
import re

##################################################
##				Configutation					##
##################################################
## write .yaml configuration filename
configfile: "/home/aquevedo/projects/test_GATKsnake/config-test.yaml"

TABLE_NAME = config["tableName"]

# Outpur directories
FASTQDIR = config["fastqdir"]
DATADIR = config["datadir"]
RESDIR = config["resdir"]
LOGDIR = config["logdir"]

# Software directories
GATK_FOLDER = config["gatk_folder"]
PICARD_FOLDER = config["picard_folder"]

# Dir with Conda ennvironment .yaml definintions and tmporary folder
CONDADIR = config["condadir"]
TMP_FOLDER = config["tmp_folder"]

## Function to create directories unless they already exist
def tryMkdir(path):
	try:
		os.makedirs(path,exist_ok=False) # raise error if exists
	except FileExistsError:
		pass

[tryMkdir(p) for p in (DATADIR, RESDIR, LOGDIR, 
					   RESDIR + "fastQC", RESDIR + "plots", RESDIR + "varaints",
					   RESDIR + "db")]


## Genome index prefixes paths
GENOME_PATH = {
	"mm9":"",
    "mm10":"" ,
    "hg19":"",
    "grch37":"",
    "hg38":"/data_genome1/References/GenomeIndices/Human/GRCh38_no_alt_plus_hs38d1/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna",
    "grch38":"/data_genome1/References/GenomeIndices/Human/GRCh38_no_alt_plus_hs38d1/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna",
    "-":""}

REF_FASTA_DICT = {
	"hg38": "/data_genome1/References/Human/Sequences/Genome/GRCh38_no_alt_plus_hs38d1_Verily.fna",
	"grch38": "/data_genome1/References/Human/Sequences/Genome/GRCh38_no_alt_plus_hs38d1_Verily.fna",
	"mm10": ""
}

DBSNP_DICT = {
	"hg38": "/data_genome1/SharedSoftware/GATK/resources_hg38/dbsnp_146.hg38.vcf",
	"grch38": "/data_genome1/SharedSoftware/GATK/resources_hg38/dbsnp_146.hg38.vcf",
	"mm10": "-"
}

GNOMAD_SITES_DICT = {
	"hg38": "/data_genome1/SharedSoftware/GATK/resources_Mutect2/af-only-gnomad.hg38.vcf.gz",
	"grch38": "/data_genome1/SharedSoftware/GATK/resources_Mutect2/af-only-gnomad.hg38.vcf.gz",
	"mm10" : "-"
}

GOLD_INDELS_DICT = {
	"hg38": "/data_genome1/SharedSoftware/GATK/resources_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
	"grch38": "/data_genome1/SharedSoftware/GATK/resources_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
	"mm10": "-"
}

PON_DICT = {
	"hg38": "-",
	"grch38": "-",
	"mm10": "-"
}

# Number chromosomes not including X (Y is discarded)
NCHR_DIC = {
	"hg38": 22,
	"hg19": 22,
	"grch38": 22,
	"grch37": 22,
	"mm10": 19,
	"mm9": 19
}

##################################################
## 				    METADATA					##
##################################################
data = pd.read_csv(TABLE_NAME,sep="\t")

def field_from_sample(Sample, field):
	# ?P<group_name> for clarity
	match = re.match(r"(?P<readID>\w+)_(?P<readLetter>\w+)_(?P<readRun>\w+)_(?P<readLane>\w+)", Sample)

	return match.group(field)

## Remove .fastq.gz to use basename with expand() in rule "all"
data["R1Basename"] = [f.replace(".fastq.gz","") for f in data["R1"]]
data["R2Basename"] = [f.replace(".fastq.gz","") for f in data["R2"]]

# Alignment
data["PathGenome"] = [GENOME_PATH[i] for i in data.Genome]

data["ReadID"] = [field_from_sample(sample, "readID") for sample in data.Samples]
data["ReadLetter"] = [field_from_sample(sample, "readLetter") for sample in data.Samples]
data["ReadRun"] = [field_from_sample(sample, "readRun") for sample in data.Samples]
data["ReadLane"] = [field_from_sample(sample, "readLane") for sample in data.Samples]

# Base quality score recalibration parame
data["RefFASTA"] = [REF_FASTA_DICT[i] for i in data.Genome]
data["GoldIndels"] = [GOLD_INDELS_DICT[i] for i in data.Genome]
data["DbSNP"] = [DBSNP_DICT[i] for i in data.Genome]

# Mutec2 resources
data["PanelNormals"] = [PON_DICT[i] for i in data.Genome]
data["NumChr"] = [NCHR_DIC[i] for i in data.Genome]

# Map Control and Tumor samples
# /* TO DO */ 
# First do it by hand, then thinnk if possible to aotomation
#	- Sample = readID_readLetter_readRun_readLane
#	- Individual = str To name mutec output
# 	- IsControl: bool
#	- NormalSample: if data.isControl == "yes"; data.ControlSample
#	- MutecOutName: if data.IsControl == "no"; data.Individual + "_somatic.vcf.gz"

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
		expand(DATADIR + "align/{sample}_rg_dedup.bam", sample = data.Samples.unique()),
		expand(DATADIR + "align/stats/remove_duplicates_{sample}.txt", sample = data.Samples.unique())


def get_readPair(pairID, fq_list):
	search = [re.match(r".+_{}_.+".format(pairID), fq) for fq in fq_list]
	filtered = list(filter(lambda x: x != None, search))

	return " ".join([filt[0] for filt in filtered])
	

rule bwa_map:
	input:
		R1 = lambda wildcards: expand(FASTQDIR + '{R1}',
			R1 = data.R1[data.Samples == wildcards.sample].values),
		R2 = lambda wildcards: expand(FASTQDIR + '{R2}',
			R2 = data.R2[data.Samples == wildcards.sample].values)
	output:
		bam = DATADIR + 'align/{sample}_sorted.bam'
	log:
		LOGDIR + 'bwa/{sample}.log'
	threads:
		get_resource("bwa", "threads") + get_resource("samtools_sort", "threads")
	resources:
		mem_mb = get_resource("bwa", "mem_mb"),
		walltime = get_resource("bwa", "walltime")
	params:
		R1 = lambda wildcards, input: get_readPair("R1", input.R1),
		R2 = lambda wildcards, input: get_readPair("R2", input.R2),
		bwa_ix = lambda wildcards: expand(data.PathGenome[data.Samples == wildcards.sample].values[0]),
		bwa_threads = get_resource("bwa", "threads"),
		sort_threads = get_resource("samtools_sort", "threads")

	shell:
		'''
		(
			bwa mem -t {params.bwa_threads} \
					   {params.bwa_ix} \
			           <(zcat {params.R1}) \
			           <(zcat {params.R2}) | \
				samtools view -h -@ 3 -F 4 -O bam - | \
				samtools sort -@ {params.sort_threads} - > {output.bam}
		) 3>&2 2>&1 1>&3 | tee -a {log}
		'''

rule add_readGroup:
	input:
		bam = DATADIR + 'align/{sample}_sorted.bam'
	output:
		rg_sorted_bam = DATADIR + 'align/{sample}_rg.bam'
	log:
		LOGDIR + 'readGroup/{sample}.log'
	threads:
		4
	resources:
		mem_mb = get_resource("samtools","mem_mb"),
		walltime = get_resource("samtools","walltime")
	params:
		RG = lambda wildcards: expand(data.RG[data.Samples == wildcards.sample]),
		PLATFORM = lambda wildcards: expand(data.PLATFORM[data.Samples == wildcards.sample])
	shell:
		'''
		(
			samtools view -F 8 -u -O BAM {input.bam} | \
			samtools addreplacerg -r \
			"@RG\tID:{params.RG}\tPL:{params.PLATFORM}"\
			 -@ 4 -O BAM -o {output.rg_sorted_bam} - 
		) 3>&2 2>&1 1>&3 | tee {log}
		'''

rule remove_duplicates:
	input:
		rg_sorted_bam = DATADIR + 'align/{sample}_rg.bam'
	output:
		nodup_bam = DATADIR + "align/{sample}_rg_dedup.bam",
		metrics = DATADIR + "align/stats/remove_duplicates_{sample}.txt"
	threads: 
		1
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime")
	params:
		tmp = TMP_FOLDER,
		gatk_folder = GATK_FOLDER
	conda:
		CONDADIR + "gatk-4.2.2.0.yaml"
	log:
		LOGDIR + "gatk/markDup_{sample}.log"
	shell:
		'''
		{params.gatk_folder}gatk \
		--java-options "-Xmx{resources.mem_mb}M" \
		--tmp-dir {params.tmp} \
		MarkDuplicates \
		--REMOVE_DUPLICATES true \
		--CREATE_INDEX true \
		-I {input.rg_sorted_bam} \
		-O {output.nodup_bam} \
		-M {output.metrics} |& tee {log}
		'''

rule createBQSR_before:
	input:
		nodup_bam = DATADIR + "align/{sample}_rg_dedup.bam"
	output:
		recal_tab = DATADIR + "align/{sample}_BSQR_before.table"
	threads: 
		1
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime")
	params:
		gatk_folder = GATK_FOLDER,
		tmp = TMP_FOLDER,
		ref_fasta = lambda wildcards: expand("{ref_fasta}",
			    ref_fasta = data.RefFASTA[data.Samples == wildcards.sample]),
		gold_indels = lambda wildcards: expand("{gold_indels}",
			    gold_indels = data.GoldIndels[data.Samples == wildcards.sample]),
		db_snp = lambda wildcards: expand("{db_snp}",
			    db_snp = data.DbSNP[data.Samples == wildcards.sample])
	conda:
		CONDADIR + "gatk-4.2.2.0.yaml"
	log:
		LOGDIR + 'gatk/createBQSR_before_{sample}.log'
	shell:
		'''
		{params.gatk_folder}gatk \
		--java-options "-Xmx{resources.mem_mb}M" \
		--tmp-dir {params.tmp} \
		BaseRecalibrator \
		-I {input.nodup_bam} \
		-R {params.ref_fasta} \
		--known_sites {params.gold_indels} \
		--known_sites {params.db_snp} \
		-O {output.recal_tab} |& tee {log}
		'''

rule applyBQSR:
	input:
		nodup_bam = DATADIR + "align/{sample}_rg_dedup.bam",
		recal_tab = DATADIR + "align/{sample}_BSQR.table"
	output:
		recal_bam = DATADIR + "align/{sample}_rg_dedup_recal.bam"
	threads:
		1
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime")
	params:
		gatk_folder = GATK_FOLDER,
		tmp = TMP_FOLDER,
		ref_fasta = lambda wildcards: expand("{ref_fasta}",
			    ref_fasta = data.RefFASTA[data.Samples == wildcards.sample])
	conda:
		CONDADIR + "gatk-4.2.2.0.yaml"
	log:
		LOGDIR + 'gatk/applyBQSR_{sample}.log'
	shell:
		'''
		{params.gatk_folder}gatk \
		--java-options "-Xmx{resources.mem_mb}M" \
		--tmp-dir {params.tmp} \
		ApplyBSQR \
		-R {params.ref_fasta} \
		-I {input.nodup_bam} \
		--bsqr-recal-file {input.recal_tab}
		-O {output.recal_bam} |& tee {log}
		'''

rule createBQSR_after:
	input:
		recal_bam = DATADIR + "align/{sample}_rg_dedup_recal.bam"
	output:
		recal_tab = DATADIR + "align/{sample}_BSQR_after.table"
	threads: 
		1
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime")
	params:
		gatk_folder = GATK_FOLDER,
		tmp = TMP_FOLDER,
		ref_fasta = lambda wildcards: expand("{ref_fasta}",
			    ref_fasta = data.RefFASTA[data.Samples == wildcards.sample]),
		gold_indels = lambda wildcards: expand("{gold_indels}",
			    gold_indels = data.GoldIndels[data.Samples == wildcards.sample]),
		db_snp = lambda wildcards: expand("{db_snp}",
			    db_snp = data.DbSNP[data.Samples == wildcards.sample])
	conda:
		CONDADIR + "gatk-4.2.2.0.yaml"
	log:
		LOGDIR + 'gatk/createBQSR_after_{sample}.log'
	shell:
		'''
		{params.gatk_folder}gatk \
		--java-options "-Xmx{resources.mem_mb}M" \
		--tmp-dir {params.tmp} \
		BaseRecalibrator \
		-I {input.nodup_bam} \
		-R {params.ref_fasta} \
		--known_sites {params.gold_indels} \
		--known_sites {params.db_snp} \
		-O {output.recal_tab} |& tee {log}
		'''


rule analyzeCovariates:
	input:
		recal_before = DATADIR + "align/{sample}_BSQR_before.table",
		recal_afer = DATADIR + "align/{sample}_BSQR_after.table"
	output:
		plots = RESDIR + "plots/{sample}_analyzeCovariates.pdf",
		csv = LOGDIR + "gatk/analyzeCovariates_csv_{sample}.log"
	threads:
		1
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime")
	params:
		gatk_folder = GATK_FOLDER,
		tmp = TMP_FOLDER
	conda:
		CONDADIR + "gatk-4.2.2.0.yaml"
	log:
		LOGDIR + "gatk/analyzeCovariates_{sample}.log"
	shell:
		'''
		{params.gatk_folder}gatk \
			--java-options "-Xmx{resources.mem_mb}M" \
			--tmp-dir {params.tmp} \
			AnalyzeCovariates
			-before {input.recal_before} \
			-after {input.recal_after} \
			-plots {output.plots} \
			-csv {output.csv} |& tee {log}
		'''




























