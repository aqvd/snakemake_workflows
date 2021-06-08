#!/bin/bash

# Extract fragment length prediction using macs2 predictd command and use it
# for calling peaks afterwards.
# macs2 built binary must be included in $PATH
#
# USAGE:
# ./predict_fragment.sh <treatBAM> <inputBAM>  <hs | mm> <sampleName> 
# <inputSampleName> <outdir> <datadir>

TREAT=$1
INPUT=$2
ORG=$3
SAMPLE=$4
INPUT_SAMPLE=$5
OUTDIR=$6
DATADIR=$7

## heck correct number of arguments
if [[ $# -ne 7 ]]
then
	echo "USAGE:  
	./predict_fragment.sh <treatBAM> <inputBAM>  <hs | mm> <sampleName> 
	<inputSampleName> <outdir> <datadir>

predict_fragment.sh requires this 7 parameters:

· BAMtreat: chip treatment

· BAMinput: chip input

· hs | mm: organism

· sampleName: to name macs2 files

· inputSampleName: name for input samples  

· outdir: if not exists, it will be created

· datadir: where alignment directory is

"
	exit 01
fi

## Check correct organism 
if [[ "${ORG}" != "hs" || "${ORG}" != "mm" ]]
then
	echo 'Second argument must be either:
	1) hs for "Homo sapiens" samples
	2) mm for "Mus musculus" samples
Using hs  as default'
	ORG=hs

fi

## Create outdir if it does not already exists
mkdir -p "${OUTDIR}"

## Predict fragment size and extract legth from stdErr redirected file
echo "===== Predicting fragment size ======="
macs2 predictd -i "${TREAT}" -g ${ORG} 2> "${OUTDIR}/${SAMPLE}_predictd.txt" &&
echo "Preciction succeded"
FragLen=$(grep -Eio "predicted fragment length is .+" \
"${OUTDIR}/${SAMPLE}_predictd.txt" | grep -Eo "[0-9]+")

echo "Fragment length is: ${FragLen}"

## Check if need for subsampling
echo "Running downsample.sh ..."

/home/aquevedo/snakemake_workflows/Chip-Seq/scripts/downsample.sh \
"${DATADIR}align/stats/picardMarkDup_${SAMPLE}.txt" \
"${DATADIR}align/stats/picardMarkDup_${INPUT_SAMPLE}.txt" \
"${DATADIR}align/stats/downsample_${SAMPLE}.txt" &&

echo "Completed downsample.sh"

## Call peaks with MACS2

DOWNSAMPLE_RES=$(grep No_need_downsampling \
${DATADIR}align/stats/downsample_${SAMPLE}.txt)

echo "${DOWNSAMPLE_RES}"

if [[ -n "${DOWNSAMPLE_RES}" ]]; then
	echo ">>No need for downsampling"

	echo "===== Calling Peaks ====="
	macs2 callpeak -t "${TREAT}" -c "${INPUT}" --gsize "${ORG}" \
	--extsize "${FragLen}" --nomodel --name "${SAMPLE}" --outdir "${OUTDIR}" &&
	
	echo "Peak calling succeded, Exiting"

	exit 00

else
	DOWN_INPUT=$(grep "input" ${DATADIR}align/stats/downsample_${SAMPLE}.txt)
	if [[ -z "${DOWN_INPUT}" ]]; then ## input not present in results
		echo ">> Downsample treatment ${SAMPLE}"
		
		FRACTION=$(sed -En '2p' ${DATADIR}align/stats/downsample_${SAMPLE}.txt)
		echo "Retained ${FRACTION} of reads"

		echo "=== DownsampleSam ==="
		gatk DownsampleSam --java-options "-Xmx4096M" \
		-I "${TREAT}" \
		-O "${TREAT/%.bam/_${SAMPLE}_downsampled.bam}" \
		-P "${FRACTION}" \
		--METRICS_FILE "${DATADIR}align/stats/DownsampleSamPicard_${SAMPLE}.txt" \
		--CREATE_INDEX true &&

		echo "Downsampling finished"
		
		echo "===== Calling Peaks ====="
		macs2 callpeak -t "${TREAT/%.bam/_${SAMPLE}_downsampled.bam}" \
		-c "${INPUT}" --gsize "${ORG}" \
		--extsize "${FragLen}" --nomodel --name "${SAMPLE}" \
		--outdir "${OUTDIR}" &&
		
		echo "Removing temporary downsampled treatment .bam:
		${TREAT/%.bam/_${SAMPLE}_downsampled.bam} "
		rm "${TREAT/%.bam/_${SAMPLE}_downsampled.bam}" "${TREAT/%.bam/_${SAMPLE}_downsampled.bai}"

		echo "Peak calling succeded, Exiting"

	else
		echo ">> Downsample input ${INPUT_SAMPLE}"
		
		FRACTION=$(sed -En '2p' ${DATADIR}align/stats/downsample_${SAMPLE}.txt)
		echo "Retained ${FRACTION} of reads"

		echo "=== DownsampleSam ==="
		gatk DownsampleSam --java-options "-Xmx4096M" \
		-I "${INPUT}" \
		-O "${INPUT/%.bam/_${SAMPLE}_downsampled.bam}" \
		-P "${FRACTION}" \
		--METRICS_FILE "${DATADIR}align/stats/DownsampleSamPicard_${INPUT_SAMPLE}_accordingTo_${SAMPLE}.txt" \
		--CREATE_INDEX true &&

		echo "Downsampling finished"
		
		echo "===== Calling Peaks ====="
		macs2 callpeak -t "${TREAT}" \
		-c "${INPUT/%.bam/_${SAMPLE}_downsampled.bam}" \
		--gsize "${ORG}" \
		--extsize "${FragLen}" --nomodel --name "${SAMPLE}" \
		--outdir "${OUTDIR}" &&

		echo "Removing temporary downsampled input .bam:
		${INPUT/%.bam/_${SAMPLE}_downsampled.bam} "
		rm "${INPUT/%.bam/_${SAMPLE}_downsampled.bam}" "${INPUT/%.bam/_${SAMPLE}_downsampled.bai}"
	
		echo "Peak calling succeded, Exiting"

	fi

fi
exit 00

