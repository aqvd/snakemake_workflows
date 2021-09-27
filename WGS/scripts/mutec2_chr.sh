#!/bin/bash

usage="Usage:
	
./mutec2_chr.sh <ref_fa>	 reference fasta
				<tumor>	 recalibrated .bam tumor
				<normal>	 racalibrated .bam normal
				<gnomad>	 .vcf gnomad sites
				<pon>	 .vcf panel of mutations in normals
				<indiv>	 individual ID to label files
				<gatk_dir> 	directory with gatk script
				<n_chroms>	 number chromosomes not including X
				<vcf_dir> 	directory for .vcf files
				<db_dir> 	directory for models, tables and intermediate files 
				<threads>	 number of jobs created by GNU parallel

Notes: ¡¡¡ Directories must end with '/' !!!
"

[[ $# -ne 8 ]]; then
	echo "${usage}"; exit 01
fi

ref_fa=$1
tumor=$2
normal=$3
gnomad=$4
pon=$5
indiv=$6
gatk_dir=$7
n_chr=$8
vcf_dir=$9
db_dir=${10}
threads=${11}

export ref_fa
export tumor
export normal
export gnomad
export pon
export indiv
export gatk_dir
export vcf_dir
export db_dir
export threads

all_chroms="$(printf "%s " $(seq 1 ${n_chr})) X"

function run_mutec2 {
	local chr=$1
	# ref_fa=$1
	# tumor=$2
	# normal=$3
	# gnomad=$4
	# pon=$5
	# indiv=$6
	# gatk_dir=$7
	# vcf_dir=$8
	# db_dir=$9

	local out_mutec2="${vcf_dir}${indiv}_chr${chr}_unfilt.vcf.gz"
	local output_f1r2="${db_dir}${indiv}_chr${chr}-f1r2.tar.gz"

	# 1- run Mutec2 with --f1r2 argument to detect strand bias later
	${gatk_dir}/gatk Mutec2 \
		-R "${ref_fa}" \
		-I "${tumor}" \
		-I "${normal}" \
		-L "${chr}" \
		--f1r2-tar-gz "${output_f1r2}" \
		-normal "${normal}" \
		--germline-resource "${gnomad}" \
		--panel-of-normals "${pon}" \
		-O "${out_mutec2}" 
}

# Function must be exported to use parallel 
export -f run_mutec2

# 2- Input all f2r2 files to model orientation bias. 
#all_fqr2_input=`for chr in ${all_chroms}; do printf -- "-I ${db_dir}${indiv}_chr${chr}-fir2.tar.gz "; done`

echo "${all_chroms}" | sed -E -e 's/ /\n/g' | parallel -j ${threads} run_mutec2 :::: ${chr}

