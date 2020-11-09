#!/bin/bash

## ./all_unique_summits.sh {outFilename} {sample1}_summits.bed {sample2}_summits.bed [{sampleN}_summits.bed]
## ./all_unique_summits.sh {outFilename} *_summits.bed  
##
## Generates a file called "summits_unique.bed" given a list of summits.bed files
## at the same directory where {sample1}_summits.bed is located.
##
## For each {sample}_summits.bed a {sample}_peaks.narrowPeak wit the same
## value of {sample} must exist. 
## Those files are the output of macs2 callpeak, so macs2 output directory should contain
## all the files necessary

USAGE="./all_unique_summits.sh {outFilename} {sample1}_summits.bed {sample2}_summits.bed [{sampleN}_summits.bed]
./all_unique_summits.sh {outFilename} *_summits.bed  "

if [ "$#" == "0" ]
then
	echo "${USAGE}"
	exit 01
fi

DIR=$(dirname $2)
RES="${DIR}/{1}"

#touch "${RES}"

## Counter to know how many files were precessed
declare -i count=1

## While loop with shift to process all files sequentially
## (( $# )) is true until there are files to be processed
while (( $# ))
do
	## replace extension in $1 to match peak file name
	PEAKS="${1/_summits.bed/_peaks.narrowPeak}"
	

	if [[ $count -eq 1 ]]
	## First time just cat sample1 summits into file called outFilename=${RES}
	then
		echo ">> Creating output file ${RES} with summits form ${2}"
		cat ${2} > ${RES}
		PEAKS="${2/_summits.bed/_peaks.narrowPeak}"

		## increase and shift
		((count+=1))
		shift 

	else if [[ $count -eq 2 ]]
	## we have already processed sample1 file in the first shift. shift to skip
		shift
	else
	## Process the rest of files
		## Save summit filename in new variable as $1 will suffer substitution
		SUMMITS="${1}"
		## replace extension in $1 to match peak file name
		PEAKS="${1/_summits.bed/_peaks.narrowPeak}"

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
