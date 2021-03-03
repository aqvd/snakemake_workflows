#!/bin/bash

TREAT=$1
INPUT=$2
ORG=$3
SAMPLE=$4
OUTDIR=$5

## heck correct number of arguments
if [[ $# -ne 5 ]]
then
	echo "USAGE:  
	./predict_fragment.sh <BAMtreat> <BAMinput> <hs | mm> <sampleName> <outdir>

predict_fragment.sh requires this 5 parameters:
BAMtreat: chip treatment
BAMinput: chip input
hs | mm: organism
sampleName: to name macs2 files
outdir: if not exists, it will be created"
	exit 01
fi

## Check correct organism 
if [[ "$ORG" != "hs" || "$ORG" != "mm" ]]
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
macs2 predictd -i "${TREAT}" -g "$ORG" 2> "${OUTDIR}/${SAMPLE}_predictd.txt" &&
echo "Preciction succeded"
FragLen=$(grep -Eio "predicted fragment length is .+" \
"${OUTDIR}/${SAMPLE}_predictd.txt" | grep -Eo "[0-9]+")

echo "Fragment length is: ${FragLen}"

## Call peaks with MACS2
echo "===== Calling Peaks ====="
macs2 callpeak -t "${TREAT}" -c "${INPUT}" --gsize "${ORG}" \
--extsize "${FragLen}" --nomodel --name "${SAMPLE}" --outdir "${OUTDIR}" &&
echo "Peak calling succeded, Exiting"
exit 00