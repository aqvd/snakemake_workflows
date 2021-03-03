for f in *.txt
do
	echo "============ ${f%%.txt} =================="
	res=$(cut -d$'\t' -f8 ${f} | gsed -E 's/[ ].+//g' | sort | uniq -c)
	tot=$(echo "${res}" | gawk 'BEGIN{s=0}; {s+=$1}; END {print s}')
	echo ">>TOTAL ANNOTATIONS: ${tot}"
	echo "${res}" | gawk -F' ' -v total=${tot} 'BEGIN{ s=0 }; \
    { p=($1/total)*100; print p,$2 }'
done