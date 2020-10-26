#!/bin/bash

# Extract fragment length prediction using macs2 predictd command
# macs2 built binary must be included in $PATH
#
# USAGE:
# ./predict_fragment.sh <BAMfile> <hs | mm> <sampleName> <outdir>

BAM=$1
ORG=$2
SAMPLE=$3
OUTDIR=$4

## heck correct number of arguments
if [[ $# -ne 4 ]]
then
	echo "USAGE: 
	./predict_fragment.sh <BAMfile> <hs | mm> <sampleName> <outdir>"
	exit 01
fi

## Check correct organism 
if [[ "$ORG" -ne "hs" ]] || [[ "$ORG" -ne "mm" ]]
then
	echo 'Second argument must be either:
	1) hs for "Homo sapiens" samples
	2) mm for "Mus musculus" samples
\	Using hs  as default'
	ORG=hs
fi

## Create outdir if it does not already exists
mkdir -p "${OUTDIR}"

macs2 predictd -i $BAM -g $ORG 2> "${OUTDIR}/${SAMPLE}_predictd.txt"

FragLen=$(grep -Eio "predicted fragment length is .+" \
"${OUTDIR}/${SAMPLE}_prediction.txt" | grep -Eo "[0-9]+")

echo "${FragLen}"
