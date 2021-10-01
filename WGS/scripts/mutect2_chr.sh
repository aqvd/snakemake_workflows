#!/bin/bash

usage="Usage:
	
./mutect2_chr.sh 
				<ref_fa>	 reference fasta

				<tumor_I>	 recalibrated .bam tumor
				<normal_I>	 racalibrated .bam normal
				<normal_sample>	 sample name (SN readgroup) nonrmal

			5	<gnomad>	 .vcf gnomad sites
				<pon>	 .vcf panel of mutations in normals
				<regions> 	ExomeSeq regions [PATH | \"chr_paralell\"]

				<indiv>	 individual ID to label files
				<n_chroms>	 number chromosomes not including X

			10	<threads>	 number of jobs created by GNU parallel
				<mem_gatk> 	mem in MB given to GATK

				<vcf_dir> 	directory for .vcf files
				<gatk_dir> 	directory with gatk script
				<db_dir> 	directory for models, tables and intermediate files 
			15	<tmp_dir> 	tmp dir for GATK

Notes: ¡¡¡ Directories must end with '/' !!!
"

if [[ $# -lt 15 ]]; then ## !!!!!!!!!!!!!!!!!!!!!!! put number arguments
	echo "Number arguments: $# < 15"
	echo "${usage}"; exit 01
fi

ref_fa=$1

tumor_I=$2
normal_I=$3
normal_sample=$4

gnomad=$5
pon=$6
regions=$7

indiv=$8
n_chr=$9

threads=${10}
mem_gatk=${11}

vcf_dir=${12}
gatk_dir=${13}
db_dir=${14}
tmp_dir=${15}

export ref_fa

export tumor_I
export normal_I
export normal_sample

export gnomad
export pon

export indiv

export threads
export mem_gatk

export gatk_dir
export vcf_dir
export db_dir

# all_chroms="$(printf "chr%s " $(seq 1 ${n_chr})) chrX"
all_chroms="$(printf "chr%s " $(seq 10 11)) chrX"

function run_mutect2 {
	local chr=$1

	echo "ref_fa = ${ref_fa}" 

	echo "tumor = ${tumor_I}" 
	echo "normal = ${normal_I}" 
	echo "normal_sample = ${normal_sample}" 

	echo "gnomad = ${gnomad}" 
	echo "pon = ${pon}" 

	echo "indiv = ${indiv}" 
	
	echo "gatk_dir = ${gatk_dir}" 
	echo "vcf_dir = ${vcf_dir}" 
	echo "db_dir = ${db_dir}"

	local out_mutect2="${vcf_dir}${indiv}_${chr}_unfilt.vcf.gz"
	local output_f1r2="${db_dir}${indiv}_${chr}-f1r2.tar.gz"

	echo "out_mutect2 = ${out_mutect2}"
	echo "output_f1r2 = ${output_f1r2}"

	# 1- run Mutect2 with --f1r2 argument to detect strand bias later
	echo -e "Starting Mutect2: Individual: ${indiv}:${chr}"

	command="${gatk_dir}/gatk Mutect2 \
		--java-options \"-Xmx${mem_gatk}M -Djava.io.tmpdir=${tmp_dir}\" \
		-R \"${ref_fa}\" \
		${tumor_I} \
		${normal_I} \
		-normal ${normal_sample} \
		-L \"${chr}\" \
		--f1r2-tar-gz \"${output_f1r2}\" \
		--germline-resource \"${gnomad}\" \
		--panel-of-normals \"${pon}\" \
		-O \"${out_mutect2}\"" #&&
	
	echo ${command}
	
	echo -e "\n Finised Mutect2"
	# echo -e " >> FilterMutectCalls Individual: ${indiv}:${chr}"

	# local out_filter_mutect2="${vcf_dir}${indiv}_${chr}_filtered.vcf.gz"
	
	# ${gatk_dir}/gatk FilterMutectCalls -R "${ref_fa}" \
	# 	-V "${out_mutect2}" \
	# 	-O "${out_filter_mutect2}" &&

	# echo -e " << Finished FilterMutectCalls"
}

# Function must be exported to use parallel 
export -f run_mutect2

if [[ "${regions}" == "chr_paralell" ]]; then
	
	echo -e "\n        -> Runing in parallel in all chromosomes: ${all_chroms}"
	echo "${all_chroms}" | sed -E -e 's/ /\n/g' | parallel -j ${threads} run_mutect2 &&
	echo -e "\n...Finised parallel mutect2..."

	# 2- Concatenate per chr results
	#    2.1- LearnReadOrientation model for all chromosomes
	all_f1r2_input=`for chr in ${all_chroms}; do 
		printf -- "-I ${db_dir}${indiv}_${chr}-f1r2.tar.gz "; done`

	echo -e "\n\t>> LearnReadOrientationModel\nINPUT FILES:\n${all_f1r2_input}"
		${gatk_dir}/gatk LearnReadOrientationModel 
			$all_f1_r2_input \
			-O ${db_dir}${indiv}_read-orientation-model.tar.gz &&
	echo -e "<< LearnReadOrientationModel Finised\n"

	#    2.2- 
	all_mutect_files=`for chr in ${all_chroms}; do 
		printf -- "${vcf_dir}${indiv}_${chr}_unfilt.vcf.gz "; done`

	echo -e "\n\t>> bcftools merge VCF files\nINPUT FILES:\n${all_f1r2_input}"
		bcftools merge ${all_mutect_files} -Oz -o "${vcf_dir}${indiv}_somatic.vcf.gz" &&
	echo -e "<< Merge VCF files Finised\n"

	exit 0

elif [[ -e ${regions} ]]; then
	
	echo -e "\n        -> Runing mutect2 in exome regions ${regions}"
	run_mutect2 ${regions} 
	echo -e "\n...Finised Mutect2"

	exit 0

fi

exit 2




