>intermediate.tmp
>all_dups_stats.tsv

for f in picardMarkDup*
do
	echo ${f##picardMarkDup}

	gawk -F "\t" -v filename=${f##picardMarkDup_} 'BEGIN {OFS="\t"; ORS=""; eol="\n"}
	/^Unknown/ { print filename,$2,$5,$6,$9,eol }' ${f} > tmp_stats.tmp

	cat all_dups_stats.tsv tmp_stats.tmp > intermediate.tmp
	cat intermediate.tmp > all_dups_stats.tsv
done

rm *.tmp
