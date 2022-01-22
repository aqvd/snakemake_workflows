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
				<f1r2_out>	name of final genome-wide f1r2 orientatino bias
				<logdir>	log directory

Notes: ¡¡¡ Directories must end with '/' !!!
"

args=("$@")

if [[ $# -lt 18 ]]; then ## !!!!!!!!!!!!!!!!!!!!!!! put number arguments
	echo "Less arguments ($#) than needed (17)"
	for ((i = 0; i < $#; i++)); do 
		echo "Arg ${i}	->	${args[$i]}"
	done
	echo "${usage}"; exit 01
fi

ref_fa=$1

tumor_I=$2
normal_I=$3
normal_sample=$4

# If not empty, add parameter name (--param) for gnomad and panel of normals
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
f1r2_out=${17}

logdir=${18}

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
export logdir

all_chroms="$(printf "chr%s " $(seq 1 ${n_chr}))chrX"
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

	local tmp_unfilt="${vcf_dir}${indiv}_${chr}_unfilt.vcf.gz"
	local tmp_f1r2="${db_dir}${indiv}_${chr}_f1r2.tar.gz"

	echo "tmp_unfilt = ${tmp_unfilt}"
    echo "tmp_f1r2 = ${tmp_f1r2}"  

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
		--f1r2-tar-gz "${tmp_f1r2}" \
		-O "${tmp_unfilt}" &&
	
	echo -e "\n Finised Mutect2"
}

# Function must be exported to use parallel 
export -f run_mutect2

if [[ "${regions}" == "chr_parallel" ]]; then
	## ============================= RUN MUTECT ===================== ##
	echo -e "\n        -> Runing in parallel in all chromosomes: ${all_chroms}\n"

	# Create also one logfile per chromosome
	echo "${all_chroms}" | sed -E -e 's/ /\n/g' \
	|parallel -j ${threads} "run_mutect2 {} > ${logdir}gatk/mutect2/${indiv}_{}.log 2>&1" &&
	
	echo -e "\n...Finised parallel mutect2..."

	# =================== Merge chr results ================== ##
	#  Merge f1r2 
	all_f1r2_input=`for chr in ${all_chroms}; do 
		printf -- "-I ${db_dir}${indiv}_${chr}_f1r2.tar.gz "; done`
	
	${gatk_dir}gatk --java-options "-Xmx${mem_gatk}M -Djava.io.tmpdir=${tmp_dir}" \
		LearnReadOrientationModel \
		$all_f1r2_input \
		-O "${f1r2_out}"

	# Merge unfiltered VCFs
	all_mutect_files=`for chr in ${all_chroms}; do 
		printf -- "-I ${vcf_dir}${indiv}_${chr}_unfilt.vcf.gz "; done`

	echo -e "\n\t>> GATK merge VCF files\nINPUT FILES:\n${all_mutect_files}"
		${gatk_dir}gatk --java-options "-Xmx${mem_gatk}M -Djava.io.tmpdir=${tmp_dir}" \
		MergeVcfs \
		${all_mutect_files} \
		-O "${vcf_out}"  &&
	echo -e "<< Merge VCF files Finised\n"

	# Merge vcf.gz.stats
	all_stats_files=`for chr in ${all_chroms}; do 
		printf -- "-stats ${vcf_dir}${indiv}_${chr}_unfilt.vcf.gz.stats "; done`
	
	echo -e "\n\t>> GATK merge VCF"
	${gatk_dir}gatk --java-options "-Xmx${mem_gatk}M -Djava.io.tmpdir=${tmp_dir}" \
		MergeMutectStats \
		${all_stats_files} \
		-O "${vcf_out}.stats"	

	#	2.4 Index (tbi) the resulting vcf
	${gatk_dir}gatk --java-options "-Xmx${mem_gatk}M -Djava.io.tmpdir=${tmp_dir}" \
		IndexFeatureFile \
		-I "${vcf_out}"
	echo -e "<< Merge VCF Finised and indexed \n"

	# remove tmp .vcf files, either without or with .stats or .tbi extensions 
	parallel rm {1}{2} ::: $all_mutect_files ::: "" ".tbi" ".stats"

	parallel rm ${db_dir}{1}_{2}_f1r2.tar.gz ::: ${indiv} ::: ${all_chroms}
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
		--f1r2-tar-gz "${f1r2_out}" \
		-O "${vcf_out}"
	
	echo -e "\n...Finised Mutect2"

	exit 0

fi

exit 2




