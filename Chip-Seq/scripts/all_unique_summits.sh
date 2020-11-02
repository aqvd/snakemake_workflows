#!/bin/bash

## ./all_unique_summits.sh 	{sample1}_summits.bed {sample2}_summits.bed [{sampleN}_summits.bed]
## ./all_unique_summits.sh 	*_summits.bed  
##
## Generates a file called "summits_unique.bed" given a list of summits.bed files
## at the same directory where {sample1}_summits.bed is located.
##
## For each {sample}_summits.bed a {sample}_peaks.narrowPeak wit the same
## value of {sample} must exist. 
## Those files are the output of macs2 callpeak, so macs2 output directory should contain
## all the files necessary

USAGE="./all_unique_summits.sh 	{sample1}_summits.bed {sample2}_summits.bed [{sampleN}_summits.bed]
./all_unique_summits.sh  *_summits.bed  "

if [ "$#" == "0" ]
then
	echo "${USAGE}"
	exit 01
fi

DIR=$(dirname $1)
RES="${DIR}/summits_unique.bed"

touch "${RES}"

## Counter to know how many files were precessed
declare -i count=1

## While loop with shift to process all files sequentially
## (( $# )) is true until there are files to be processed
while (( $# ))
do
	## Save summit filename in new variable as $1 will suffer substitution
	SUMMITS="${1}"
	## replace extension in $1 to match peak file name
	PEAKS="${1/_summits.bed/_peaks.narrowPeak}"
	

	if [[ $count -eq 1 ]]
	then
		echo ">> Processing file ${count} ${SUMMITS}"
		PEAKS2="${2/_summits.bed/_peaks.narrowPeak}"
		
		echo -e "Peaks Files are: \n ${PEAKS} \n ${PEAKS2}"
		
		## intersect first 2 files, get summits for peaks unique in -a
		bedtools intersect -v -a "${PEAKS}" -b "${PEAKS2}" > int1.tmp
		echo "REGIONS BEFORE: $(wc -l int1.tmp | sed -E 's/[ a-zA-Z.].+//g')"
		
		bedtools intersect -wa -a "${SUMMITS}" -b int1.tmp >> "${RES}"
		echo "Intersected regions: $(wc -l ${SUMMITS} | sed -E 's/[ a-zA-Z.].+//g')"
		echo "REGIONS AFTER $(wc -l ${RES} | sed -E 's/[ a-zA-Z.].+//g')"
		
		## increase and shift
		((count+=1))
		shift 

	else
		echo ">>file ${count} ${SUMMITS}"
		echo -e "Peaks file is: \n ${PEAKS}"
		
		## The same logoc as in the if clause
		## Save peaks not already present to int1.tmp and get their summits
		bedtools intersect -v -a "${PEAKS}" -b "${RES}" > int1.tmp &&
		bedtools intersect -wa -a "${SUMMITS}" -b int1.tmp > int2.tmp &&

		
		## count lines befor and after appending new summits to ${RES}
		echo "REGIONS BEFORE: $(wc -l ${RES} | sed -E 's/[ a-zA-Z.].+//g')"
		cat int2.tmp >> "${RES}"
		
		echo "Intersected regions $(wc -l ${SUMMITS} | sed -E 's/[ a-zA-Z.].+//g')"
		echo "REGIONS AFTER $(wc -l ${RES} | sed -E 's/[ a-zA-Z.].+//g')"
		## increase counter and shift
		((count+=1))
		shift 

	fi
done

echo "...Remonvig temporary files"
rm int1.tmp int2.tmp

echo ">>> FINISHED. Exiting... <<<"
exit 00
