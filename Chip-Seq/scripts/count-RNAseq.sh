#!/usr/bin/bash

## count-RNAseq.sh <gtf> <type_feature> <ID_countMatrix> <is_stranded> 
##				   <output_countMatrix> <bam.1> [<bam.2> ... <bam.n>]
##
## Script to call htseq-count in Snakefile and extract countMAtrix column
## names from <bam.x> filenames.
## IMPORTANT! .bam files MUST BE the last arguments

gtf=$1
type_feature=$2
ID_countMatrix=$3
is_stranded=$4
output_countMatrix=$5
bams=( ${@:5:${#}} ) # arguments from 5 onwards

USAGE=""
if [[ "${#}" -le 6  ]]; then
	echo "USAGE:
	count-RNAseq.sh <gtf> <type_feature> <ID_countMatrix> <is_stranded> 
				   <output_countMatrix> <bam.1> [<bam.2> ... <bam.n>]

Not enough arguments: ${#} 
Please specify the 5 mandatory first ones and at least one .bam file"
	exit 01
fi

## extract bam basenames for naming columns. 1st column gene_ids, so whitespace
echo "Provided alignment .bam files:
${bams[@]}"

colnames=(" ")

for bam in ${bams[@]}; do
	sample=${bam##/*/}
	colnames+=( "${sample/.bam/}" )
done

echo -e "Colnames are:\n${colnames[@]}"

## Generate 1st line of countMatrix with column names. Replace spaces by tabs. 
echo "${colnames[@]}" | gsed -E 's/ /\t/g' > "${output_countMatrix}"

## Run htseq-count and append its output to 1st line
htseq_command="htseq-count --format=bam --order=pos --stranded=${is_stranded} \
--type=${type_feature} --idattr=${ID_in_countMatrix} \
--mode=union --secondary-alignments=ignore \
${bams[@]} ${gtf} >> ${output_countMatrix}"

echo -e "Running htseq-count command...\n${htseq_command}"
${htseq_command} &&

echo "htseq-count finished. 
Exit..."

exit 00


