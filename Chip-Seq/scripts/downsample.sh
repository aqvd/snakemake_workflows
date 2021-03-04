dupStats=$1
input_dupStats=$2

out_scale=$3

validReads=$(awk -F "\t" '/^Unknown/ {total=$2; dup=$6} 
END {print total-dup}' ${dupStats})

echo "Valid reads treatment: ${validReads}"

input_validReads=$(awk -F "\t" '/^Unknown/ {total=$2; dup=$6} 
END {print total-dup}' ${input_dupStats})

echo "Valid reads input: ${input_validReads}"

scale=$(python -c "print( float(${validReads}) / float(${input_validReads}) )")

echo "Scale Factor: treatment / input = ${scale}"

if (( $(echo "${scale} > 0.95" | bc -l) & $(echo "${scale} < 1.05" | bc -l) )); then
	echo "No_need_downsampling" > ${out_scale} 
	exit 0
fi

if (( $(echo "${scale} < 0.95" | bc -l) )); then # more reads in input
	echo "Downsample_Control: ${input_dupStats##/*/}" > ${out_scale}
	echo "scale factor for downsampling is: ${scale}"

	result=$(python -c "print( round(${scale},2) )")
	echo ${result} >> ${out_scale} 
	exit 0

else # more reads tratment
	echo "Downsample_Treatmemt: ${dupStats##/*/}" > ${out_scale}
	scale=$(python -c "print( 1/${scale} )")
	echo "scale factor for downsampling is: ${scale}"

	result=$(python -c "print( round(${scale},2) )")
	echo ${result} >> ${out_scale} 
	exit 0
fi


