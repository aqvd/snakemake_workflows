## 
## USAGE: 
## bwMean.sh R1.bw R2.bw outfile
## 
USAGE="USAGE: 
	bwMean.sh R1.bw R2.bw outfile"

chromSizes="/home/aquevedo/opt/homer/data/genomes/hg19/chrom.sizes"
#chromSizes="/Users/aqo/Desktop/ALVARO/Software/homer/data/genomes/hg19/chrom.sizes"
WD=$(dirname ${3})

if [[ $# -ne 3 ]]
then
	echo "${USAGE}"
	exit 01
fi

echo "Obtaining mean score of ${1##.+/} and ${2##.+/}"
wiggletools mean ${1} ${2} > ${WD}/mean.tmp

echo "Sorting intermediary .bg file "
bedSort ${WD}/mean.tmp ${WD}/mean.tmp &&

echo "From BedGraph to BigWig"
bedGraphToBigWig ${WD}/mean.tmp ${chromSizes} ${3} &&

echo "${3##./} bigWig generated at ${WD}"
ls -lh ${3}

echo "Removing temporary file"
rm ${WD}/mean.tmp

echo "Finished. Exiting..."

exit 00