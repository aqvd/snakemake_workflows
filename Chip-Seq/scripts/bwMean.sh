## 
## USAGE: 
## bwMean.sh R1.bw R2.bw outfile
## 
USAGE="USAGE: 
	bwMean.sh R1.bw R2.bw outfile"

chromSizes="/home/aquevedo/opt/homer/data/genomes/hg19/chrom.sizes"
#chromSizes="/Users/aqo/Desktop/ALVARO/Software/homer/data/genomes/hg19/chrom.sizes"

WD=$(dirname ${3})
tempFile="${WD}/${1##.+/}_mean.tmp"

## When no need to merge replicates, outfile is a symlink
if [[ $# -eq 2 ]]
then
	echo "Only one bigWig, no need to merge. 
${2##.+/} is just a symbolic link to ${1##.+/}"
	ln -sf "${1}" "${2}"
	exit 00
fi

if [[ $# -ne 3 ]]
then
	echo "${USAGE}"
	exit 01
fi

echo "Obtaining mean score of ${1##.+/} and ${2##.+/}"
wiggletools mean ${1} ${2} > ${tempFile}

echo "Sorting intermediary .bg file "
bedSort ${tempFile} ${tempFile} &&

echo "From BedGraph to BigWig"
bedGraphToBigWig ${tempFile} ${chromSizes} ${3} &&

echo "${3##./} bigWig generated at ${WD}"
ls -lh ${3}

echo "Removing temporary file"
rm ${tempFile}

echo "Finished. Exiting..."

exit 00