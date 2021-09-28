#!/bin/bash

# When using --split-3 option of fast(er)q-dump, files are named:
# SRA-Accesion followed vy  _1.fastq _2.fastq or .fastq if no mate present
# 
# For calrity, rename it to know what is each file.

if [[ $# -ne 5 ]]; then
	echo -e "Wrong number of arguments. 5 needed, $# provided
	USAGE:
./name_split3 DIR SRA PROT COND REP"
fi

DIR=$1
SRA=$2

PROT=$3
COND=$4
REP=$5

for f in ${DIR}${SRA}*.fastq ; do
	if [ -e "${f}" ]; then
		mate_id="${f##*${SRA}_}"
		echo -e "Processing mate: ${mate_id} \n"

		new_filename="${DIR}${PROT}_${COND}_${REP}_${SRA}_R${mate_id}"

		mv "${f}" "${new_filename}" &&
		echo -e "renaming ${f##*/} to ${new_filename##*/} \n"

		gzip "${new_filename}" && 
		echo -e "compressing ${f##*/} using gzip \n"
		
	else
		echo -e "File ${f} does not exist"
		continue
	fi
done