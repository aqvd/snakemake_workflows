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
			--java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={params.tmp}" \
			Mutec2 \
			-R {params.ref_fasta} \
			-I {input.bam} \
			--max-mnp-distance 0 \
			--af-of-alleles-not-in-resource \
			-O {output} |& tee {log}
		'''

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

def repeat_argument(argument, value_list, **kwags):
	"""
	sometimes a argument (eg: -I <file>) must be specified multiple times,
	this function creates the appropriate concatenation
	"""

	argument_l = [str(argument) + str(v) for v in value_list]
	res = " ".join(argument_l)

	return(res)

def expand_argument(path, string_format,
					type_sample = ["tumor", "sample"],
					indiv, df,
				    argument, value_list, **kwags):
	"""
	generates repeated argument value pairs (eg -I path/file_1.txt -I path/file_2.txt)
	using expand function from snakemake to create the paths.

	string_format must be '{expansion}_bla_bla_bla.xxx'
	"""
	if not path.endswith("/"):
		path = path + "/"

	if type_sample == "tumor":
		arg_value_l = expand(path + string_to_expand, 
						  expansion = tumor_sample_from_individual(indiv, df))
	else:
		arg_value_l = expand(path + string_to_expand, 
						  expansion = nnormal_sample_from_individual(indiv, df))

	res = repeat_argument(argument, arg_value_l)

	return(res)



print(tumor_sample_from_individual("mouseA", data))

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
		RESDIR + "variants/{indiv}_somatic.vcf.gz",
		RESDIR + "db/{indiv}_read-orientation-model.tar.gz"
	threads:
		1
	resources:
		mem_mb = get_resource("gatk", "mem_mb"),
		walltime = get_resource("gatk","walltime"),
	params:
		reference = lambda wildcards: data.IxPrefPath[data.Individual == wildcards.indiv].values[0],

		tumor_I_param = lambda wildcards: 
			expand_argument(DATADIR + "align", "{expansion}_rg_dedup_recal.bam",
						    data, "Samples", "", "-I ", Individual = "indivA", IsControl = "no"),
		normal_I_param = lambda wildcards: 
			expand_argument(DATADIR + "align", "{expansion}_rg_dedup_recal.bam",
						    data, "Samples", "", "-I ", Individual = "indivA", IsControl = "yes"),
		normal_name = lambda wildcards:
			get_column_df(data, "Samples", "", Individual = wildcards.indiv, IsControl = "yes"),

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
		{params.scrpit_folder}mutect2_chr.sh \
			{params.reference} \
			{params.tumor_I_param} \
			{params.normal_I_param} \
			{params.normal_name} \
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