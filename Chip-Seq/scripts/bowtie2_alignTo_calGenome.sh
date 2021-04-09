#!/bin/bash

## ========= bowtie2_alignTo_calGenome.sh ========
## 
## Script that executes shell commands for rule "bowtie2_alignTo_calGenome"
## so that conda: directive is allowed (not alowed if run directive)
##
## All the arguments but the reads are unique strings. To deal with variable
## number of reads, pass them in the last place and use array slice (0 based)
## to access al of them.

args=( "${@}" )

out_bam=$1
out_stats=$2
stat_dir=$3

tmp_unal=$4

calGenIx=$5
genomeIndex=$6

threads=$7
log=$8
## Reads are the last arguments
reads=${args[@]:8}

if [[ "${calGenIx}" == "NoCalibration" ]]; then
	echo "No calibration"

	echo -e "reads are:\n${reads[@]}"
	echo -e "============================================\n"
	echo "Aligning against reference genome"

	(
		bowtie2 -U ${reads[@]} -x ${genomeIndex} -p ${threads} \
			--time --no-unal | \
		samtools view -h -@ 3 -O bam - > ${out_bam} 
	) 3>&2 2>&1 1>&3 | \
	tee ${log} &&

	## Create the rest of output files, but empty, to avoid missingOutputException
	mkdir -p ${stat_dir} && 
	echo "No needs calibration" > ${out_stats}
else
	echo "Yes calibration"

	echo -e "reads are:\n${reads[@]}"
	echo -e "============================================\n"
	echo "Aligning against reference genome"
	
	## Get only reads that align to reference genome: {out_bam}
	## Get reads that do NOT align to rederence: {params.tmp_unal}.
	(
		bowtie2 -x ${genomeIndex} -U ${reads[@]} -p ${threads} --time \
			--un-gz ${tmp_unal} --no-unal | \
		samtools view -h -@ 3 -O bam - > ${out_bam} 
	) 3>&2 2>&1 1>&3 | \
	tee ${log} &&

	echo "Aligning against calubration genome genome"
	## {output.stats}: stats alignments reads unique to calibration genome
	bowtie2 -x ${calGenIx} -U ${tmp_unal} \
		-p ${threads} --time --no-unal -S /dev/null |& \
	tee ${out_stats} && \
	rm ${tmp_unal}
fi

exit 00