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

				<vcf_out>	name of final genome-wide unfiltered vcf

Notes: ¡¡¡ Directories must end with '/' !!!
"

if [[ $# -lt 16 ]]; then ## !!!!!!!!!!!!!!!!!!!!!!! put number arguments
	echo "Number arguments: $# < 16"
	echo "${usage}"; exit 01
fi

ref_fa=$1

tumor_I=$2
normal_I=$3
normal_sample=$4

# If not empty parameters for gnomad and panel of normals
gnomad=${5:+"--germline-resource ${5}"}
pon=${6:+"--panel-of-normals ${6}"}

regions=$7

indiv=$8
n_chr=$9

threads=${10}
mem_gatk=${11}

vcf_dir=${12}
gatk_dir=${13}
db_dir=${14}
tmp_dir=${15}

vcf_out=${16}

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
export tmp_dir

all_chroms="$(printf "chr%s " $(seq 1 ${n_chr})) chrX"
# all_chroms="$(printf "chr%s " $(seq 10 11)) chrX"

function run_mutect2 {
	local chr=$1

	echo "ref_fa = ${ref_fa}" 
	echo "regions = ${chr}"

	echo "tumor = ${tumor_I}" 
	echo "normal = ${normal_I}" 
	echo "normal_sample = ${normal_sample}" 

	echo "${gnomad}" 
	echo "${pon}" 

	echo "indiv = ${indiv}" 
	
	echo "gatk_dir = ${gatk_dir}" 
	echo "vcf_dir = ${vcf_dir}" 
	echo "tmp_dir = ${tmp_dir}" 
	echo "db_dir = ${db_dir}"

	local out_unfilt="${vcf_dir}${indiv}_${chr}_unfilt.vcf.gz"
# !	local output_f1r2="${db_dir}${indiv}_${chr}-f1r2.tar.gz" !

	echo "out_unfilt = ${out_unfilt}"
# !  echo "output_f1r2 = ${output_f1r2}"  

	# 1- run Mutect2 with --f1r2 argument to detect strand bias later
	echo -e "\n\t >> ======= Starting Mutect2: Individual: ${indiv}:${chr} ========= <<"

	${gatk_dir}gatk --java-options "-Xmx${mem_gatk}M -Djava.io.tmpdir=${tmp_dir}" \
		Mutect2 \
		-R "${ref_fa}" \
		${tumor_I} \
		${normal_I} \
		${normal_sample} \
		-L "${chr}" \
		${gnomad} \
		${pon} \
		-O "${out_unfilt}" &&
	
	echo -e "\n Finised Mutect2"
}

# Function must be exported to use parallel 
export -f run_mutect2

if [[ "${regions}" == "chr_parallel" ]]; then
	
	echo -e "\n        -> Runing in parallel in all chromosomes: ${all_chroms}\n"
	echo "${all_chroms}" | sed -E -e 's/ /\n/g' | parallel -j ${threads} run_mutect2 &&
	echo -e "\n...Finised parallel mutect2..."

	# 2- Concatenate per chr results

	#    2.2- Merge unfiltered VCFs
	all_mutect_files=`for chr in ${all_chroms}; do 
		printf -- "${vcf_dir}${indiv}_${chr}_unfilt.vcf.gz "; done`

	echo -e "\n\t>> bcftools merge VCF files\nINPUT FILES:\n${all_mutect_files}"
		bcftools merge ${all_mutect_files} -O z -o "${vcf_out}"  &&
	echo -e "<< Merge VCF files Finised\n"

	exit 0

elif [[ -e ${regions} ]]; then

	echo "ref_fa = ${ref_fa}" 
	echo "regions = ${regions}"

	echo "tumor = ${tumor_I}" 
	echo "normal = ${normal_I}" 
	echo "normal_sample = ${normal_sample}" 

	echo "${gnomad}" 
	echo "${pon}" 

	echo "indiv = ${indiv}" 
	
	echo "gatk_dir = ${gatk_dir}" 
	echo "vcf_dir = ${vcf_dir}" 
	echo "tmp_dir = ${tmp_dir}" 
	echo "db_dir = ${db_dir}"
	
	echo -e "\n        -> Runing mutect2 in exome regions ${regions}\n"
	
	command="${gatk_dir}gatk \
		--java-options \"-Xmx${mem_gatk}M -Djava.io.tmpdir=${tmp_dir}\" \
		Mutect2 \
		-R \"${ref_fa}\" \
		${tumor_I} \
		${normal_I} \
		${normal_sample} \
		-L \"${regions}\" \
		${gnomad} \
		${pon} \
		-O \"${vcf_out}\""

	echo $command

	${gatk_dir}gatk --java-options "-Xmx${mem_gatk}M -Djava.io.tmpdir=${tmp_dir}" \
		Mutect2 \
		-R "${ref_fa}" \
		${tumor_I} \
		${normal_I} \
		${normal_sample} \
		-L "${regions}" \
		${gnomad} \
		${pon} \
		-O "${vcf_out}"
	
	echo -e "\n...Finised Mutect2"

	exit 0

fi

exit 2




