import os
import pandas as pd
import numpy as np
import glob
import re
from utils_snakemake_workflows.utils import *

##################################################
##				Configutation					##
##################################################
## write .yaml configuration filename
configfile: "/home/aquevedo/snakemake_workflows/test_WGS/config-test.yaml"

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
SCRIPT_FOLDER = config["script_folder"]

## Function to create directories unless they already exist
def tryMkdir(path):
	try:
		os.makedirs(path,exist_ok=False) # raise error if exists
	except FileExistsError:
		pass

[tryMkdir(p) for p in (DATADIR + "stats", LOGDIR, TMP_FOLDER,
					   RESDIR + "fastQC", RESDIR + "plots", RESDIR + "variants",
					   RESDIR + "db")]


## Genome index prefixes paths
GENOME_IX_PREFIX_DICT = {
	"mm9":"/home/aquevedo/genomes/mouse/mm38/GRCm38.p6_standardAndMito_bwa",
	"mm38": "/home/aquevedo/genomes/mouse/mm38/GRCm38.p6_standardAndMito_bwa",
    "mm10":"" ,
    "mm39": "",
    "hg19":"",
    "grch37":"",
    "hg38":"/home/aquevedo/resources/references/GRCh38/Verily_decoy/bwa-0.7.12/GRCh38_Verily_v1.genome.fa",
    "grch38":"/home/aquevedo/resources/references/GRCh38/Verily_decoy/bwa-0.7.12/GRCh38_Verily_v1.genome.fa",
    "-":""}

REF_FASTA_DICT = {
	"mm9":"/home/aquevedo/genomes/mouse/mm38/GRCm38.p6_standardAndMito.fa",
	"mm38": "/home/aquevedo/genomes/mouse/mm38/GRCm38.p6_standardAndMito.fa",
    "mm10":"/home/aquevedo/genomes/mouse/mm39/mm39.fa" ,
    "mm39": "/home/aquevedo/genomes/mouse/mm39/mm39.fa",
    "hg19":"",
    "grch37":"",
	"hg38": "/home/aquevedo/resources/references/GRCh38/Verily_decoy/GRCh38_Verily_v1.genome.fa",
	"grch38": "/home/aquevedo/resources/references/GRCh38/Verily_decoy/GRCh38_Verily_v1.genome.fa",
}

DBSNP_DICT = {
	"mm9":"/home/aquevedo/genomes/mouse/mm38/dbSNP_UCSC-chr-names.vcf.gz",
	"mm38": "/home/aquevedo/genomes/mouse/mm38/dbSNP_UCSC-chr-names.vcf.gz",
    "mm10":"" ,
    "mm39": "",
	"hg38": "/data_genome1/SharedSoftware/GATK/resources_hg38/dbsnp_146.hg38.vcf",
	"grch38": "/data_genome1/SharedSoftware/GATK/resources_hg38/dbsnp_146.hg38.vcf",
}

GNOMAD_SITES_DICT = {
	"mm9":"/home/aquevedo/genomes/mouse/mm38/mgp_REL2005_snps_f1-8_AF_chrUCSC_afOnly_bialellelic.sorted.vcf.gz",
	"mm38": "/home/aquevedo/genomes/mouse/mm38/mgp_REL2005_snps_f1-8_AF_chrUCSC_afOnly_bialellelic.sorted.vcf.gz",
    "mm10":"" ,
    "mm39": "",
	"hg38": "/data_genome1/SharedSoftware/GATK/resources_Mutect2/af-only-gnomad.hg38.vcf.gz",
	"grch38": "/data_genome1/SharedSoftware/GATK/resources_Mutect2/af-only-gnomad.hg38.vcf.gz",
}

GOLD_INDELS_DICT = {
	"mm9":"/home/aquevedo/genomes/mouse/mm38/mgp_REL2005_snps_f1-8_AF_chrUCSC_afOnly_bialellelic.sorted.vcf.gz",
	"mm38": "/home/aquevedo/genomes/mouse/mm38/mgp_REL2005_snps_f1-8_AF_chrUCSC_afOnly_bialellelic.sorted.vcf.gz",
    "mm10":"" ,
    "mm39": "",
	"hg38": "/data_genome1/SharedSoftware/GATK/resources_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
	"grch38": "/data_genome1/SharedSoftware/GATK/resources_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
}

PON_DICT = {
	"mm9":"",
	"mm38": "",
    "mm10":"" ,
    "mm39": "",
	"hg38": "/home/aquevedo/resources/GATK_public_data/1000g_pon.hg38.vcf.gz",
	"grch38": "/home/aquevedo/resources/GATK_public_data/1000g_pon.hg38.vcf.gz",
}

REGIONS_MUTECT_DIC = {
	"mm9":"/home/aquevedo/genomes/mouse/mm38/S0276129_Regions.bed",
	"mm38": "/home/aquevedo/genomes/mouse/mm38/S0276129_Regions.bed",
    "mm10":"" ,
    "mm39": "",
	"hg38": "/data_genome1/References/AgilentSureSelect/Human_Exome_V6_UTR/hg38/S07604624_Covered_fixed.bed",
	"grch38": "/data_genome1/References/AgilentSureSelect/Human_Exome_V6_UTR/hg38/S07604624_Covered_fixed.bed",
	"chr_parallel": "chr_parallel"
}

EXOME_TARGETS_DICT = {
	"mm9":"/home/aquevedo/genomes/mouse/mm38/S0276129_Regions.bed",
	"mm38": "/home/aquevedo/genomes/mouse/mm38/S0276129_Regions.bed",
    "mm10":"" ,
    "mm39": "",
	"hg38": "",
	"grch38": "",
}

# Number chromosomes not including X (Y is discarded)
NCHR_DIC = {
	"hg38": 22,
	"hg19": 22,
	"grch38": 22,
	"grch37": 22,
	"mm9":19,
	"mm38": 19,
    "mm10":19 ,
    "mm39": 19,
}

##################################################
## 				    METADATA					##
##################################################
data = pd.read_csv(TABLE_NAME,sep="\t")

## Remove .fastq.gz to use basename with expand() in rule "all"
data["R1Basename"] = [f.replace(".fastq.gz","") for f in data["R1"]]
data["R2Basename"] = [f.replace(".fastq.gz","") for f in data["R2"]]

# Alignment
data["IxPrefPath"] = [GENOME_IX_PREFIX_DICT[i] for i in data.Genome]

# To add reag group in SAM files
data["ReadID"] = [field_from_sample(sample, "readID") for sample in data.SamplesIDalign]
data["ReadLetter"] = [field_from_sample(sample, "readLetter") for sample in data.SamplesIDalign]
data["ReadRun"] = [field_from_sample(sample, "readRun") for sample in data.SamplesIDalign]
data["ReadLane"] = [field_from_sample(sample, "readLane") for sample in data.SamplesIDalign]
data["RunLane"] = data.ReadRun + "_" + data.ReadLane
data["ReadGroupID"] = data.RunLane + "_" + data.SamplesIDalign

# Base quality score recalibration parame
data["RefFASTA"] = [REF_FASTA_DICT[i] for i in data.Genome]
data["GoldIndels"] = [GOLD_INDELS_DICT[i] for i in data.Genome]
data["Gnomad"] = [GNOMAD_SITES_DICT[i] for i in data.Genome]
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
		# expand(DATADIR + "align/{sample}_unsorted.bam", sample = data.SamplesIDalign.unique())
		# expand(DATADIR + "align/{sample}_rg_dedup.bam", sample = data.SamplesIDalign.unique()),
		# expand(DATADIR + 'align/{sample}_rg_dedup_merged.bam',sample = data.Samples.unique())
		# expand(DATADIR + "align/{sample}_BSQR_before.table", sample = data.Samples.unique()),
		# expand(DATADIR + "align/{sample}_BSQR_after.table", sample = data.Samples.unique()),
		# expand(RESDIR + "plots/{sample}_analyzeCovariates.pdf", sample = data.Samples.unique())
		# expand(RESDIR + "variants/{indiv}_somatic_unfilt.vcf.gz", indiv = data.Individual.unique()),
		# expand(RESDIR + "db/{indiv}_read-orientation-model.tar.gz",indiv = data.Individual.unique()),
		expand(RESDIR + "variants/{indiv}_somatic_filtered.vcf.gz", indiv = data.Individual.unique()),


def get_readPair(pairID, fq_list):
	"""
	In case several .fastq per read mate, get all to concatenate with zcat prior bwa
	"""
	search = [re.match(r".+_{}_.+".format(pairID), fq) for fq in fq_list]
	filtered = list(filter(lambda x: x != None, search))

	return " ".join([filt[0] for filt in filtered])
	

rule bwa_map:
	input:
		R1 = lambda wildcards: expand(FASTQDIR + '{R1}',
			R1 = data.R1[data.SamplesIDalign == wildcards.sample].values),
		R2 = lambda wildcards: expand(FASTQDIR + '{R2}',
			R2 = data.R2[data.SamplesIDalign == wildcards.sample].values)
	output:
		bam = temp(DATADIR + 'align/{sample}_unsorted.bam')
	log:
		LOGDIR + 'bwa/{sample}.log'
	threads:
		get_resource("bwa", "threads") + 5
	resources:
		mem_mb = get_resource("bwa", "mem_mb"),
		walltime = get_resource("bwa", "walltime")
	params:
		R1 = lambda wildcards, input: get_readPair("R1", input.R1),
		R2 = lambda wildcards, input: get_readPair("R2", input.R2),
		bwa_ix = lambda wildcards: data.IxPrefPath[data.SamplesIDalign == wildcards.sample].values[0],
		bwa_threads = get_resource("bwa", "threads") - 5
	shell:
		'''
		(
			bwa mem -t {params.bwa_threads} \
				{params.bwa_ix} \
				<(zcat {params.R1}) \
				<(zcat {params.R2}) | \
			samtools view -h -@ 5 -O BAM - > {output.bam}
		) 3>&2 2>&1 1>&3 | tee -a {log}
		'''

rule sort_bam:
	input:
		DATADIR + 'align/{sample}_unsorted.bam'
	output:
		bam = DATADIR + 'align/{sample}_sorted.bam'
	log:
		LOGDIR + 'gatk/SortSam_{sample}.log'
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime")
	params:
		tmp = TMP_FOLDER,
		gatk_folder = GATK_FOLDER
	shell:
		'''
		{params.gatk_folder}gatk \
		--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
			SortSam \
			-I {input} \
			-O {output} \
			--SORT_ORDER coordinate \
			-VALIDATION_STRINGENCY LENIENT |& tee {log} && 
       
        {params.gatk_folder}gatk \
		--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
		BuildBamIndex \
			-I {output} \
			-VALIDATION_STRINGENCY LENIENT
		'''


rule add_readGroup:
	input:
		bam = DATADIR + 'align/{sample}_sorted.bam'
	output:
		rg_sorted_bam = DATADIR + 'align/{sample}_rg.bam'
	log:
		LOGDIR + 'samtools/readGroup_{sample}.log'
	threads:
		4
	resources:
		mem_mb = get_resource("samtools","mem_mb"),
		walltime = get_resource("samtools","walltime")
	params:
		ID = lambda wildcards: expand(data.ReadGroupID[data.SamplesIDalign == wildcards.sample]),
		PL = lambda wildcards: expand(data.Platform[data.SamplesIDalign == wildcards.sample]),
		PU = lambda wildcards: expand(data.RunLane[data.SamplesIDalign == wildcards.sample]),
		LB = lambda wildcards: wildcards.sample,
		SM = lambda wildcards: wildcards.sample
	shell:
		'''
		# -F 12 to filter unmaped and mate unmaped
		# view GATK AddReplaceReadGRoup to understand the tags
		(
			samtools view -@ 3 -F 12 -u -O BAM {input.bam} | \
			samtools addreplacerg -r \
				"@RG\tID:{params.ID}\tPL:{params.PL}\tPU:{params.PU}\tLB:{params.LB}\tSM:{params.SM}"\
				-@ {threads} -O BAM - > {output.rg_sorted_bam}
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
		--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
		MarkDuplicates \
		--REMOVE_DUPLICATES true \
		--CREATE_INDEX true \
		-I {input.rg_sorted_bam} \
		-O {output.nodup_bam} \
		-M {output.metrics} |& tee {log}
		'''

rule mergeBAM:
	input:
		lambda wildcards: expand(DATADIR + "align/{sampIDalign}_rg_dedup.bam",
			sampIDalign = data.SamplesIDalign[wildcards.sample == data.Samples])
	output:
		mergedBAM = DATADIR + 'align/{sample}_rg_dedup_merged.bam'
	threads:
		1
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime")
	params:
		tmp = TMP_FOLDER,
		gatk_folder = GATK_FOLDER,
		ip_str = lambda wildcards, input: " -I ".join(input)
	log:
		LOGDIR + 'gatk/merge_{sample}.log'
	shell:
		'''
		{params.gatk_folder}gatk \
		--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
		MergeSamFiles \
			-I {params.ip_str} \
			-O {output.mergedBAM} 
		'''


rule createBQSR_before:
	input:
		merged_bam = DATADIR + "align/{sample}_rg_dedup_merged.bam"
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
			    ref_fasta = data.RefFASTA[data.Samples == wildcards.sample])[0],
		gold_indels = lambda wildcards: expand("{gold_indels}",
			    gold_indels = data.GoldIndels[data.Samples == wildcards.sample])[0],
		db_snp = lambda wildcards: expand("{db_snp}",
			    db_snp = data.DbSNP[data.Samples == wildcards.sample])[0]
	conda:
		CONDADIR + "gatk-4.2.2.0.yaml"
	log:
		LOGDIR + 'gatk/createBQSR_before_{sample}.log'
	shell:
		''' 
		#
		if [ -n "{params.gold_indels}" ]; then
			gold_indels="--known-sites "{params.gold_indels}
		else
			gold_indels=" "
		fi

		if [ -n "{params.db_snp}" ]; then
			db_snp="--known-sites "{params.db_snp}
		else
			db_snp=" "
		fi

		{params.gatk_folder}gatk \
		--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
		BaseRecalibrator \
		-I {input.merged_bam} \
		-R {params.ref_fasta} \
		${{gold_indels}} \
		${{db_snp}} \
		-O {output.recal_tab} |& tee {log}
		'''

rule applyBQSR:
	input:
		nodup_bam = DATADIR + "align/{sample}_rg_dedup_merged.bam",
		recal_tab = DATADIR + "align/{sample}_BSQR_before.table"
	output:
		recal_bam = DATADIR + "align/{sample}_rg_dedup_merged_recal.bam"
	threads:
		1
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime")
	params:
		gatk_folder = GATK_FOLDER,
		tmp = TMP_FOLDER,
		ref_fasta = lambda wildcards: expand("{ref_fasta}",
			    ref_fasta = data.RefFASTA[data.Samples == wildcards.sample])[0]
	conda:
		CONDADIR + "gatk-4.2.2.0.yaml"
	log:
		LOGDIR + 'gatk/applyBQSR_{sample}.log'
	shell:
		'''
		{params.gatk_folder}gatk \
		--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
		ApplyBQSR \
		-R {params.ref_fasta} \
		-I {input.nodup_bam} \
		--bqsr-recal-file {input.recal_tab} \
		-O {output.recal_bam} |& tee {log}
		'''

rule createBQSR_after:
	input:
		recal_bam = DATADIR + "align/{sample}_rg_dedup_merged_recal.bam"
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
			    ref_fasta = data.RefFASTA[data.Samples == wildcards.sample])[0],
		gold_indels = lambda wildcards: expand("{gold_indels}",
			    gold_indels = data.GoldIndels[data.Samples == wildcards.sample])[0],
		db_snp = lambda wildcards: expand("{db_snp}",
			    db_snp = data.DbSNP[data.Samples == wildcards.sample])[0]
	conda:
		CONDADIR + "gatk-4.2.2.0.yaml"
	log:
		LOGDIR + 'gatk/createBQSR_after_{sample}.log'
	shell:
		'''
		if [ -n "{params.gold_indels}" ]; then
			gold_indels="--known-sites "{params.gold_indels}
		else
			gold_indels=" "
		fi

		if [ -n "{params.db_snp}" ]; then
			db_snp="--known-sites "{params.db_snp}
		else
			db_snp=" "
		fi

		{params.gatk_folder}gatk \
		--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
		BaseRecalibrator \
		-I {input.recal_bam} \
		-R {params.ref_fasta} \
		${{gold_indels}} \
		${{db_snp}} \
		-O {output.recal_tab} |& tee {log}
		'''

rule analyzeCovariates:
	input:
		recal_before = DATADIR + "align/{sample}_BSQR_before.table",
		recal_after = DATADIR + "align/{sample}_BSQR_after.table"
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
			--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
			AnalyzeCovariates \
			-before {input.recal_before} \
			-after {input.recal_after} \
			-plots {output.plots} \
			-csv {output.csv} |& tee {log}
		'''

rule mutec2_tumor_vs_normal:
	input:
		tumor = lambda wildcards: 
			expand(DATADIR + "align/{sample}_rg_dedup_merged_recal.bam",
				   sample = get_column_df(data, 
				   						  column = "Samples",
				   						  filt_out = "", 
				   					      Individual = wildcards.indiv,
				   					      IsControl = "no"))[0],
		normal = lambda wildcards: 
			expand(DATADIR + "align/{sample}_rg_dedup_merged_recal.bam",
				   sample = get_column_df(data, 
				   						  column = "Samples",
				   						  filt_out = "", 
				   					      Individual = wildcards.indiv,
				   					      IsControl = "yes"))[0],
	output:
		vcf = RESDIR + "variants/{indiv}_somatic_unfilt.vcf.gz",
		f1r2 = RESDIR + "db/{indiv}_f1r2.tar.gz",
		stats = RESDIR + "variants/{indiv}_somatic_unfilt.vcf.gz.stats"
	threads:
		get_resource("Mutect2", "threads")
	resources:
		mem_mb = get_resource("Mutect2", "mem_mb"),
		walltime = get_resource("Mutect2","walltime")
	params:
		reference = lambda wildcards: data.RefFASTA[data.Individual == wildcards.indiv].values[0],

		tumor_I_param = lambda wildcards: 
			expand_argument(path = DATADIR + "align/", 
							string_to_expand = "{expansion}_rg_dedup_merged_recal.bam",
						    df = data, 
						    column = "Samples", 
						    filt_out = "", 
						    argument = "-I ", 
						    Individual = wildcards.indiv, 
						    IsControl = "no"),

		normal_I_param = lambda wildcards: 
			expand_argument(path = DATADIR + "align/", 
							string_to_expand = "{expansion}_rg_dedup_merged_recal.bam",
						    df = data, 
						    column = "Samples", 
						    filt_out = "", 
						    argument = "-I ",
						    Individual = wildcards.indiv, 
						    IsControl = "yes"),

		normal_name = lambda wildcards:
			repeat_argument( "-normal ",
							get_column_df(df = data, 
									  column = "SamplesIDalign", 
									  filt_out = "", 
									  Individual = wildcards.indiv, 
									  IsControl = "yes")),

		gnomad = lambda wildcards: data.Gnomad[data.Individual == wildcards.indiv].values[0],
		pon = lambda wildcards: data.PanelNormals[data.Individual == wildcards.indiv].values[0],
		regions = lambda wildcards: data.RegionsMutect[data.Individual == wildcards.indiv].values[0],

		n_chroms = lambda wildcards: data.NumChr[data.Individual == wildcards.indiv].values[0],
		
		vcf_dir = RESDIR + "variants/",
		gatk_folder = GATK_FOLDER,
		db_dir = RESDIR + "db/",
		
		script_folder = SCRIPT_FOLDER,
		tmp = TMP_FOLDER,
		logdir = LOGDIR
	conda:
		CONDADIR + "gatk-4.2.2.0.yaml"
	log:
		LOGDIR + "gatk/mutect2/mutec2_tum_vs_norm_{indiv}.log"
	shell:
		'''
		{params.script_folder}mutect2_chr.sh \
			{params.reference} \
			"{params.tumor_I_param}" \
			"{params.normal_I_param}" \
			"{params.normal_name}" \
			"{params.gnomad}" \
			"{params.pon}" \
			{params.regions} \
			{wildcards.indiv} \
			{params.n_chroms} \
			{threads} \
			{resources.mem_mb} \
			{params.vcf_dir} \
			{params.gatk_folder} \
			{params.db_dir} \
			{params.tmp} \
			{output.vcf} \
			{output.f1r2} \
			{params.logdir} |& tee {log}
		'''

rule get_pileup_summaries_tumor:
	input:
		tumor_bam = lambda wildcards: 
		expand(DATADIR + "align/{sample}_rg_dedup_merged_recal.bam",
			   sample = data.Samples[(data.Individual == wildcards.indiv) & (data.IsControl == "no") ])
	output:
		RESDIR + "db/{indiv}_pileups.table"
	threads:
		1
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime")
	params:
		input = lambda wildcards, input: repeat_argument( "-I ", set(input)),
		gnomad = lambda wildcards: data.Gnomad[data.Individual == wildcards.indiv].values[0],
		gatk_folder = GATK_FOLDER,
		tmp = TMP_FOLDER
	log:
		LOGDIR + "gatk/pileup_summaries_tumor{indiv}.log"
	shell:
		'''
		if [ -n "{params.gnomad}" ]; then
			V_param="-V {params.gnomad}"
			L_param="-L {params.gnomad}"
		else
			V_param=""
			L_param=""
		fi

		{params.gatk_folder}/gatk \
			--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
			GetPileupSummaries \
			{params.input} \
			${{L_param}} \
			${{V_param}} \
			-O {output} |& tee {log}
		'''

rule calculate_contamination_tumor:
	input:
		RESDIR + "db/{indiv}_pileups.table"
	output:
		RESDIR + "db/{indiv}_contamination.table"
	threads:
		1
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime")
	params:
		gatk_folder = GATK_FOLDER,
		tmp = TMP_FOLDER
	log:
		LOGDIR + "gatk/calculate_contamination_tumor{indiv}.log"
	shell:
		'''
		{params.gatk_folder}/gatk \
			--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
			CalculateContamination \
			-I {input} \
			-O {output} |& tee {log}
		'''

rule filter_mutect_calls:
	input:
		vcf = RESDIR + "variants/{indiv}_somatic_unfilt.vcf.gz",
		read_orient = RESDIR + "db/{indiv}_f1r2.tar.gz",
		contam = RESDIR + "db/{indiv}_contamination.table",
		in_stats = RESDIR + "variants/{indiv}_somatic_unfilt.vcf.gz.stats"
	output:
		filt = RESDIR + "variants/{indiv}_somatic_filtered.vcf.gz",
		filt_stats = RESDIR + "variants/{indiv}_somatic_filtering_stats.vcf.gz"
	threads:
		1
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime")
	params:
		gatk_folder = GATK_FOLDER,
		tmp = TMP_FOLDER,
		reference = lambda wildcards: data.RefFASTA[data.Individual == wildcards.indiv].values[0],
	log:
		LOGDIR + "gatk/filter_mutect_calls_{indiv}.log"
	shell:
		'''
		{params.gatk_folder}gatk \
			--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
			FilterMutectCalls \
			-R {params.reference} \
			-V {input.vcf} \
			--stats {input.in_stats} \
	        --contamination-table {input.contam} \
	        --ob-priors {input.read_orient} \
	        --filtering-stats {output.filt_stats} \
	        -O {output.filt} |& tee {log}
		'''























