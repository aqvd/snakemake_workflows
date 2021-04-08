#!/usr/bin/env python3

import re
import sys
import pandas as pd

files = sys.argv[1:]

results_df=pd.DataFrame(columns=["sample", "validReads_treatment","validReads_input",
	"whichDownsampled", "scale_downsampling", "keptReads_downsamp",
	"macs_treatmentReads","macs_inputReads", "Errors"])

for file in files:
	sample=re.sub(".log","",file)

	with open(file,"r") as f:
		for line in f.readlines():
			line=line.rstrip()
			#		
			if re.search("Valid reads treatment:", line):
				validTreat=line.split(" ")[3]
				print(validTreat)
			#
			if re.search("Valid reads input:", line):
				validInput=line.split(" ")[3]
				print(validInput)
			#
			if re.search("scale factor for downsampling is:", line):
				scale=line.split(" ")[5]
				print(scale)
			#
			if re.search(">> Downsample (input|treatment) (.*)", line):
				whichDownsampled=re.sub(">> Downsample (input|treatment) (.*)",
					r'\2',line)
				print(f'whichDownsampled = {whichDownsampled}')

			if re.search("Kept [0-9]+ out of [0-9]+ reads", line):
				kept=re.sub(".*Kept ([0-9]+) out of [0-9]+ reads.*", r'\1', line)
				print(f'kept = {kept}')
			
			if re.search("tags after filtering in treatment:", line):
				macsTreat=re.sub(".*tags after filtering in treatment: ([0-9]+)",
					r'\1' , line)
				print(f'macsTreat = {macsTreat}')

			if re.search("tags after filtering in control:", line):
				macsCont=re.sub(".*tags after filtering in control: ([0-9]+)",
					r'\1', line)
				print(f'macsCont = {macsCont}')

			if re.search("error", line, flags=re.IGNORECASE):
				errors="Yes"
			else:
				errors="No"

		results_df=results_df.append({"sample": sample, 
				"validReads_treatment": validTreat, "validReads_input": validInput,
		 		"whichDownsampled":whichDownsampled, "scale_downsampling": scale,
		 		 "keptReads_downsamp": kept, "macs_treatmentReads": macsTreat, 
		 		 "macs_inputReads": macsCont, "Errors": errors}, ignore_index=True)

		f.close()

print(results_df)

results_df.to_csv('macsLog_results.tsv', sep='\t', index = False)

exit()



