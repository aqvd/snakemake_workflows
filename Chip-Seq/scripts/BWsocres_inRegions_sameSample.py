import pyBigWig as BW
import os
import glob
import numpy as np
import pandas as pd
import re

import matplotlib.pyplot as plt
import seaborn as sns
%pylab


ctcfDir = "/Users/aqo/Desktop/cluster/MCF10A/siNipbl/Peaks/merged_replicates_callpeak/SAonly-common_allPeaks/Cohesin_intersected_CTCF"
ctcf_regions = [os.path.join(ctcfDir, reg) for reg in 
	   ('cohesin_ctcf_siNipblCHIP.bed',"cohesin_no_ctcf_siNipblCHIP.bed")]

cohesinDir = "/Users/aqo/Desktop/cluster/MCF10A/siNipbl/Peaks/merged_replicates_callpeak/"
cohesin_regions = [os.path.join(cohesinDir, reg) for reg in 
	   ('positions_ctcf',"positions_no_ctcf")]

beds = ctcf_regions + cohesin_regions

peaks_siNIPBL_dir = "/Users/aqo/Desktop/cluster/MCF10A/siNipbl-siSA_combined/log2ratio_siNipbl_siSA_ChIPs_newPseudocount/sorted_regions_log2desc"
beds = [ reg for reg in glob.glob(peaks_siNIPBL_dir + "/*siNipbl*.bed")]

# wd = "/Users/aqo/Desktop/cluster/HAP1_cells_MAU2ko_WAPLko/log2Ratio_MAU2ko_and_WAPLko_vs_WT/cohesinCTCF_noCTCF"
wd = "/Users/aqo/Desktop/cluster/MCF10A/siNipbl-siSA_combined/log2ratio_siNipbl_siSA_ChIPs_newPseudocount/sorted_regions_log2desc"

os.chdir(wd)

bws = [os.path.join(os.getcwd(), reg) for reg in glob.glob("*.bw")]

# Find middle point if regions are not summits

def findCenter(start, end):
    if end - start != 1:
        start = int((end + start) / 2)
        end = start + 1 

    return start, end

# Read scores from bw using pyBigWig and cohesinCTCF noCTCF regions

def scoresFromBed(bed_file, opened_bw, upstream, downstream):
    with open(bed_file, "r") as bed:
        chroms=[]
        starts=[]
        ends=[]
        scores=[]
        for line in bed:
        	if line.startswith("#"):
        		continue

           	regions=line.rstrip().split("\t")
            try:
	            scores.append(opened_bw.stats(regions[0], int(regions[1]),
	                                          int(regions[2]), type="mean"))
	        except RuntimeError:
	        	continue

            chroms.append(regions[0])
            start, end = findCenter(int(regions[1]), int(regions[2]))
            starts.append(start - upstream)
            ends.append(end + upstream)

        bed.close()

        return(np.array(chroms),
               np.array(starts, dtype=np.uint),
               np.array(ends, dtype=np.uint),
               np.ma.masked_array(scores, mask=np.isnan(scores), fill_value=0))


def tmpDF(bed_file_path, ch, st, end, scores):
    df_tmp=pd.DataFrame(
        {
            "chr": ch,
            "start": st,
            "end": end,
            "regions": [os.path.basename(bed_file_path) for
               i in range(len(ch))],
            "score": scores.data.flatten()
        },

    )

    return df_tmp

def protCond_fromFile (filePath):
    file = os.path.basename(filePath)
    m = re.search(r"^(.+?)_(.+?)_.*$", file)
    prot, cond = m.group(1,2)

    return prot, cond

def check_sameSample(filePath1, filePath2):
    return protCond_fromFile(filePath1) == protCond_fromFile(filePath2)

del data
errors = 0
processed_bws = []
processed_regions = []

for bw_filename in bws:
    processed_bws.append(bw_filename)
    opened_bw = BW.open(bw_filename, "r")
    
    for bed_file in beds:
        if check_sameSample(bw_filename, bed_file):
        	print(f"\nObtaining scores from:\n{os.path.basename(bw_filename)}\n\
                At regions:\n{os.path.basename(bed_file)}")
            processed_regions.append(bed_file)
            print(f"Processed BEDs : {len(processed_regions)}")
            try:
    	        chroms, starts, ends, scores = scoresFromBed(bed_file, opened_bw, 50, 50)
            except RuntimeError:
            	errors += 1

            if len(processed_regions) == 1:
                data = tmpDF(bed_file, chroms, starts,
                               ends, scores)
                continue

            elif len(processed_regions) > 1 :
                df_tmp = tmpDF(bed_file, chroms, starts,
                                ends, scores)
                data = pd.concat([data, df_tmp], join="outer", axis=0)

        else:
            continue
    opened_bw.close()

data.shape
data.columns
data.regions.unique()
# Check number of lines in bed matches numer of rows of data
# !wc -l *siNipbl_log2Ratio*

# Missing Data as zero (np.ma.masked_array)
np.sum(data.iloc[:, 2:].values == 0)

# write to table
data.to_csv("data_log2scores_siNIPBL_sameSample.tsv", sep="\t", header=True, index= False)

