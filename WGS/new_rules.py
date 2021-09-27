def normal_sample_from_individual(indiv, df):
	df_indiv = df[df.Individual == indiv]
	res = df_indiv.loc[df_indiv["IsControl" == "yes"], "NormalSample"].values[0]

	return([res])

print(normal_sample_from_individual("mouseA", data))

rule mutec2_tumor_only:
	input:
		bam = lambda wildcards: expand(
			DATADIR + "align/{sample}_rg_dedup_recal.bam",
				sample = data.Samples[data.Individual == wildcards.indiv].values[0]
			)
	output:
		RESDIR + "variants/{indiv}_tumor_only.vcf.gz"
	threads:
		1
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime")
	params:
		gatk_folder = GATK_FOLDER,
		tmp = TMP_FOLDER,
		ref_fasta = lambda wildcards: expand("{ref_fasta}",
			ref_fasta = data.RefFASTA[data.Individual == wildcards.indiv].values[0])
	conda:
		CONDADIR + "gatk-4.2.2.0.yaml"
	log:
		LOGDIR + "gatk/mutec2_tumor_only_{indiv}.log"
	shell:
		'''
		{params.gatk_folder}gatk \
			--java-options "-Xmx{resources.mem_mb}M" \
			--tmp-dir {params.tmp} \
			Mutec2 \
			-R {params.ref_fasta} \
			-I {input.bam} \
			--max-mnp-distance 0 \
			--af-of-alleles-not-in-resource \
			-O {output} |& tee {log}
		'''

def repeat_parameter(param, value_list, **kwags):
	"""
	sometimes a parameter (eg: -I <file>) must be specified multiple times,
	this function creates the appropriate concatenation
	"""

	param_l = [str(param) + str(v) for v in value_list]
	res = " ".join(param_l)

	return(res)


## I've fund a pon file in the googleCloud from GATK. not the optimal,
## but better than nothing for my data

# rule create_ponDB:
# 	input:
# 		vcfs = expand(RESDIR + "variants/{indiv}_tumor_only.vcf.gz",
# 			indiv = data.Individual[data.IsControl == "yes"].unique())
# 	output:
# 		RESDIR + "db/pon_db_created.bnk"
# 	threads:
# 		1
# 	resources:
# 		mem_mb = get_resource("gatk", "mem_mb"),
# 		walltime = get_resource("gatk","walltime")
# 	params:
# 		gatk_folder = GATK_FOLDER,
# 		tmp = TMP_FOLDER,
# 		variants = lambda widcards, input: repeat_parameter("-V", input.vcfs)
# 		db_dir = RESDIR + "db/",
# 	conda:
# 		CONDADIR + "gatk-4.2.2.0.yaml"
# 	log:
# 		LOGDIR + "gatk/mutec2_tumor_only_{indiv}.log"
# 	shell:



def tumor_sample_from_individual(indiv, df):
	df_indiv = df[df.Individual == indiv]
	res = df_indiv.loc[df_indiv["IsControl" == "no"], "TumorSample"].values[0]

	return([res])

print(tumor_sample_from_individual("mouseA", data))

rule mutec2_tumor_vs_normal:
	input:
		tumor = lambda wildcards: expand(
			DATADIR + "align/{sample}_rg_dedup_recal.bam",
				sample = tumor_sample_from_individual(wildcards.indiv, data)
			),
		normal = lambda wildcards: expand(
			DATADIR + "align/{sample}_rg_dedup_recal.bam",
				sample = normal_sample_from_individual(wildcards.indiv, data)
			)
	output:
		RESDIR + "variants/{indiv}_somatic.vcf.gz"
	threads:
		1
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime"),
		normal_name = lambda wildcards: expand(
			"{sample}_rg_dedup_recal",
				sample = normal_sample_from_individual(wildcards.indiv, data)
			),
		pon = lambda wildcards: data.PanelNormals[data.Individual == wildcards.indiv].values[0]
	params:
		gatk_folder = GATK_FOLDER,
		tmp = TMP_FOLDER,
	conda:
		CONDADIR + "gatk-4.2.2.0.yaml"
	log:
		LOGDIR + "gatk/mutec2_tum_vs_norm_{indiv}.log"
	shell:
		'''
	
		'''