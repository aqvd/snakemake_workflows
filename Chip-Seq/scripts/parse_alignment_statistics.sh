>intermediate.tmp
>all_alignments_stats.tsv

for f in bowtie2_calibration*
do
	echo ${f##bowtie2_calibration_}

	gawk -v filename=${f##bowtie2_calibration_} 'BEGIN {eof="\n"; ORS="\t"; OFS="\t"; print filename}
	/reads; of these:/ { print $1}
	/aligned exactly 1 time/ { print $1,$2}
	/aligned >1 times/ { print $1,$2}
	/aligned 0 times/ { print $1,$2}
	END {print eof}' ${f} > tmp_stats.tmp

	cat all_alignments_stats.tsv tmp_stats.tmp > intermediate.tmp
	cat intermediate.tmp > all_alignments_stats.tsv
done

rm *.tmp