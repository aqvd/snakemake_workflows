#!/bin/bash

# Script to run in parallel Mutect2 in tumor-only mode to learn read 
# orientation bias


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