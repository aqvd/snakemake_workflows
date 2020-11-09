BEGIN{ 
	FS="\t"
	OFS="\t"
	prevTAD=0
	borderDist=0
	lastTADend=0
	lastTADchr=0
}
FNR == 1 {
	## Process first line
	prevTAD=$3
	$3=$2 + 1
	print $0
}
FNR > 1 && FNR < TNR {
	## From second to last line
	borderDist=($2-prevTAD) / 2
	$2=prevTAD + borderDist
	prevTAD=$3  
	$3=$2 + 1
	print $0
}
FNR == TNR {
	## When we arrive to last line: Save chromosome and end of TAD
	lastTADchr=$1
	lastTADend=$3
	## Calculate border as before and print it
	borderDist=($2-prevTAD) / 2
	$2=prevTAD + borderDist
	prevTAD=$3  
	$3=$2 + 1
	print $0
	## Print newline with the end of last tad
	print lastTADchr , lastTADend, lastTADend + 1
}
	