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
	"mm9":"",
    "mm10":"" ,
    "hg19":"",
    "grch37":"",
    "hg38":"/home/aquevedo/resources/references/GRCh38/Verily_decoy/bwa-0.7.12/GRCh38_Verily_v1.genome.fa",
    "grch38":"/home/aquevedo/resources/references/GRCh38/Verily_decoy/bwa-0.7.12/GRCh38_Verily_v1.genome.fa",
    "-":""}

REF_FASTA_DICT = {
	"hg38": "/home/aquevedo/resources/references/GRCh38/Verily_decoy/GRCh38_Verily_v1.genome.fa",
	"grch38": "/home/aquevedo/resources/references/GRCh38/Verily_decoy/GRCh38_Verily_v1.genome.fa",
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
	"hg38": "/home/aquevedo/resources/GATK_public_data/1000g_pon.hg38.vcf.gz",
	"grch38": "/home/aquevedo/resources/GATK_public_data/1000g_pon.hg38.vcf.gz",
	"mm10": "-"
}

REGIONS_MUTECT_DIC = {
	"hg38": "/data_genome1/References/AgilentSureSelect/Human_Exome_V6_UTR/hg38/S07604624_Covered_fixed.bed",
	"grch38": "/data_genome1/References/AgilentSureSelect/Human_Exome_V6_UTR/hg38/S07604624_Covered_fixed.bed",
	"mm10": "-"
}

EXOME_TARGETS_DICT = {
	"hg38":
	"grch38":
	"mm10":
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
data["IxPrefPath"] = [GENOME_IX_PREFIX_DICT[i] for i in data.Genome]

data["ReadID"] = [field_from_sample(sample, "readID") for sample in data.Samples]
data["ReadLetter"] = [field_from_sample(sample, "readLetter") for sample in data.Samples]
data["ReadRun"] = [field_from_sample(sample, "readRun") for sample in data.Samples]
data["ReadLane"] = [field_from_sample(sample, "readLane") for sample in data.Samples]

# Base quality score recalibration parame
data["RefFASTA"] = [REF_FASTA_DICT[i] for i in data.Genome]
data["GoldIndels"] = [GOLD_INDELS_DICT[i] for i in data.Genome]
data["Gnomad"] = [GNOMAD_SITES_DICT[i] for i in data.Genome]
data["DbSNP"] = [DBSNP_DICT[i] for i in data.Genome]

# Mutec2 resources
data["PanelNormals"] = [PON_DICT[i] for i in data.Genome]
data["NumChr"] = [NCHR_DIC[i] for i in data.Genome]
data["RegionsMutect"] = [REGIONS_MUTECT_DIC[i] for i in data.Genome]

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
		# expand(DATADIR + "align/{sample}_BSQR_before.table", sample = data.Samples.unique()),
		# expand(DATADIR + "align/{sample}_BSQR_after.table", sample = data.Samples.unique()),
		# expand(RESDIR + "plots/{sample}_analyzeCovariates.pdf", sample = data.Samples.unique())
		expand(RESDIR + "variants/{indiv}_somatic.vcf.gz", indiv = data.Individual.unique()),
		expand(RESDIR + "db/{indiv}_read-orientation-model.tar.gz",indiv = data.Individual.unique())


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
			R1 = data.R1[data.Samples == wildcards.sample].values),
		R2 = lambda wildcards: expand(FASTQDIR + '{R2}',
			R2 = data.R2[data.Samples == wildcards.sample].values)
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
		bwa_ix = lambda wildcards: data.IxPrefPath[data.Samples == wildcards.sample].values[0],
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
		LOGDIR + 'samtools/sort_{sample}.log'
	threads:
		get_resource("samtools_sort", "threads")
	resources:
		mem_mb = get_resource("samtools_sort", "mem_mb"),
		walltime = get_resource("samtools_sort", "walltime")
	shell:
		'''
		samtools sort -@ {threads} -O BAM {input} > {output.bam} 2> {log}
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
		RG = lambda wildcards: expand(data.ReadGroup[data.Samples == wildcards.sample]),
		PL = lambda wildcards: expand(data.Platform[data.Samples == wildcards.sample]),
		PU = lambda wildcards: expand(data.PlatformUnit[data.Samples == wildcards.sample]),
		LB = lambda wildcards: expand(data.Library[data.Samples == wildcards.sample]),
		SM = lambda wildcards: wildcards.sample
	shell:
		'''
		# -F 12 to filter unmaped and mate unmaped
		# view GATK AddReplaceReadGRoup to understand the tags
		(
			samtools view -@ 3 -F 12 -u -O BAM {input.bam} | \
			samtools addreplacerg -r \
				"@RG\tID:{params.RG}\tPL:{params.PL}\tPU:{params.PU}\tLB:{params.LB}\tSM:{params.SM}"\
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
		--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
		BaseRecalibrator \
		-I {input.nodup_bam} \
		-R {params.ref_fasta} \
		--known-sites {params.gold_indels} \
		--known-sites {params.db_snp} \
		-O {output.recal_tab} |& tee {log}
		'''

rule applyBQSR:
	input:
		nodup_bam = DATADIR + "align/{sample}_rg_dedup.bam",
		recal_tab = DATADIR + "align/{sample}_BSQR_before.table"
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
		--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
		ApplyBQSR \
		-R {params.ref_fasta} \
		-I {input.nodup_bam} \
		--bqsr-recal-file {input.recal_tab} \
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
		--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
		BaseRecalibrator \
		-I {input.recal_bam} \
		-R {params.ref_fasta} \
		--known-sites {params.gold_indels} \
		--known-sites {params.db_snp} \
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
			expand(DATADIR + "align/{sample}_rg_dedup_recal.bam",
				   sample = get_column_df(data, "Samples", "", 
				   					      Individual = wildcards.indiv,
				   					      IsControl = "no")),
		normal = lambda wildcards: 
			expand(DATADIR + "align/{sample}_rg_dedup_recal.bam",
				   sample = get_column_df(data, "Samples", "", 
				   					      Individual = wildcards.indiv,
				   					      IsControl = "yes")),
	output:
		RESDIR + "variants/{indiv}_{chr}_somatic_unfilt.vcf.gz",
		RESDIR + "db/{indiv}_{chr}-f1r2.tar.gz"
	threads:
		get_resource("Mutect2", "threads")
	resources:
		mem_mb = get_resource("Mutect2", "mem_mb"),
		walltime = get_resource("Mutect2","walltime")
	params:
		reference = lambda wildcards: data.RefFASTA[data.Individual == wildcards.indiv].values[0],

		tumor_I_param = lambda wildcards: 
			expand_argument(DATADIR + "align", "{expansion}_rg_dedup_recal.bam",
						    data, "Samples", "", "-I ", Individual = "indivA", IsControl = "no"),
		normal_I_param = lambda wildcards: 
			expand_argument(DATADIR + "align", "{expansion}_rg_dedup_recal.bam",
						    data, "Samples", "", "-I ", Individual = "indivA", IsControl = "yes"),
		normal_name = lambda wildcards:
			get_column_df(data, "Samples", "", Individual = wildcards.indiv, IsControl = "yes")[0],

		gnomad = lambda wildcards: data.Gnomad[data.Individual == wildcards.indiv].values[0],
		pon = lambda wildcards: data.PanelNormals[data.Individual == wildcards.indiv].values[0],
		regions = "chr_paralell",

		n_chroms = lambda wildcards: data.NumChr[data.Individual == wildcards.indiv].values[0],
		
		vcf_dir = RESDIR + "variants/",
		gatk_folder = GATK_FOLDER,
		db_dir = RESDIR + "db/",
		
		script_folder = SCRIPT_FOLDER,
		tmp = TMP_FOLDER
	conda:
		CONDADIR + "gatk-4.2.2.0.yaml"
	log:
		LOGDIR + "gatk/mutec2_tum_vs_norm_{indiv}.log"
	shell:
		'''
		{params.script_folder}mutect2_chr.sh \
			{params.reference} \
			"{params.tumor_I_param}" \
			"{params.normal_I_param}" \
			"{params.normal_name}" \
			{params.gnomad} \
			{params.pon} \
			{params.regions} \
			{wildcards.indiv}
			{params.n_chroms} \
			{threads} \
			{resources.mem_mb} \
			{params.vcf_dir} \
			{params.gatk_folder} \
			{params.db_dir} \
			{params.tmp} 	
		'''

rule learn_read_orientation_model:
	input:
		lambda wildcards: expand(RESDIR + "db/" + "{indiv}_{chr}-f1r2.tar.gz", 
			indiv = wildcards.indiv, 
			chr = ["chr" + str(x) for x in (range(1, data.NumChr[data.Individual == wildcards.indiv] + 1), "X")])
	output:
		RESDIR + "db/" + "{indiv}_read-orientation-model.tar.gz"
	threads:
		1
	resources:
		mem_mb = get_resource("Mutect2", "mem_mb"),
		walltime = get_resource("Mutect2","walltime")
	params:
		I = lambda wildcards, innput: repeat_argument("-I ", input),
		gatk_folder = GATK_FOLDER,
		tmp = TMP_FOLDER
	shell:
		'''
		{params.gatk_folder}/gatk \
			--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
			LearnReadOrientationModel \
			{params.I} \
			-O {output}		
		'''

rule get_pileup_summaries_tumor:
	input:
		tumor_bam = lambda wildcards: 
		expand(DATADIR + "align/{sample}_rg_dedup_recal.bam",
			   sample = data.Samples[(data.Individual == wildcads.indiv) & (data.IsControl == "no") ])
	output:
		RESDIR + "db/{indiv}_pileups.table"
	threads:
		1
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime")
	params:
		gnomad = lambda wildcards: data.Gnomad[data.Individual == wildcards.indiv].values[0],
		gatk_folder = GATK_FOLDER,
		tmp = TMP_FOLDER
	shell:
		'''
		{params.gatk_folder}/gatk \
			--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
			GetPileupSummaries \
			-I {input.tumor_bam} \
			-L {params.gnomad} \
			-V {params.gnomad} \
			-O {output}
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
	shell:
		'''
		{params.gatk_folder}/gatk \
			--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
			CalculateContamination \
			-I {input} \
			-O {output}
		'''

rule filter_mutect_calls:
	input:
		vcf = RESDIR + "variants/{indiv}_{chr}_somatic_unfilt.vcf.gz",
		read_orient = RESDIR + "db/{indiv}_{chr}-f1r2.tar.gz",
		contam = RESDIR + "db/{indiv}_contamination.table"
	output:
		RESDIR + "variants/{indiv}_{chr}_somatic_filtered.vcf.gz"
	threads:
		1
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime")
	params:
		gatk_folder = GATK_FOLDER,
		tmp = TMP_FOLDER
	shell:
		'''
		{params.gatk_folder}/gatk \
			--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
			FilterMutectCalls \
			-V {input.vcf} \
	        --contamination-table {input.contam} \
	        --ob-priors {input.read_orient} \
	        -O {output}
		'''























