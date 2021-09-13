import pyBigWig as BW
import os
import glob
import numpy as np
import pandas as pd
import re

import matplotlib.pyplot as plt
import seaborn as sns
%pylab


cohesinDir = "/Users/aqo/Desktop/cluster/MCF10A/siNipbl/Peaks/merged_replicates_callpeak/SAonly-common_allPeaks/Cohesin_intersected_CTCF"
cohesin_regions = [os.path.join(cohesinDir, reg) for reg in 
	   ('cohesin_ctcf_siNipblCHIP.bed',"cohesin_no_ctcf_siNipblCHIP.bed")]

ctcfDir = "/Users/aqo/Desktop/cluster/MCF10A/siNipbl/Peaks/merged_replicates_callpeak/"
ctcf_regions = [os.path.join(ctcfDir, reg) for reg in 
	   ('positions_ctcf',"positions_no_ctcf")]

# beds = ctcf_regions + cohesin_regions

peaks_siNIPBL_dir = "/Users/aqo/Desktop/cluster/MCF10A/siNipbl-siSA_combined/log2ratio_siNipbl_siSA_ChIPs_newPseudocount/sorted_regions_log2desc"
beds = [ reg for reg in glob.glob(peaks_siNIPBL_dir + "/*siNipbl*.bed")]

beds = ctcf_regions
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


def tmpDF(bed_file_path, ch, st, end, scores, score_columnName):
    df_tmp=pd.DataFrame(
        {
            "chr": ch,
            "start": st,
            "end": end,
            "regions": [os.path.basename(bed_file_path) for
               i in range(len(ch))],
            score_columnName: scores.data.flatten()
        },

        # index=[r + '_' + str(i) for r,i in zip(
        # 	[ os.path.basename(bed_file_path) for i in range(len(ch)) ], 
        # 	range(len(ch)))],
    )

    return df_tmp


del data
processed_bws = []
for bw_filename in bws:
    processed_bws.append(bw_filename)
    opened_bw = BW.open(bw_filename, "r")
    score_columnName = os.path.basename(bw_filename) + "_score"
    processed_regions = []
    errors = 0
    for bed_file in beds:
    	print(f"\nObtaining scores from:\n{os.path.basename(bw_filename)}\nAt regions:\n{os.path.basename(bed_file)}")
        processed_regions.append(bed_file)
        try:
	        chroms, starts, ends, scores = scoresFromBed(bed_file, opened_bw, 50, 50)
        except RuntimeError:
        	errors += 1

        if len(processed_regions) == 1:
            df_tmp = tmpDF(bed_file, chroms, starts,
                           ends, scores, score_columnName)
            continue

        elif len(processed_regions) > 1 and len(processed_regions) < len(beds):
            df_tmp2 = tmpDF(bed_file, chroms, starts,
                            ends, scores, score_columnName)
            df_tmp = pd.concat([df_tmp, df_tmp2], join="outer", axis=0)
            continue

        elif len(processed_regions) == len(beds):
            df_tmp2 = tmpDF(bed_file, chroms, starts,
                            ends, scores, score_columnName)
            df_tmp = pd.concat([df_tmp, df_tmp2], join="outer", axis=0)

            if len(processed_bws) == 1:
            	print("{:^80}".format("-- bigWig 1, no merge"))
                data = df_tmp
                data.index = pd.Index(range(data.shape[0]))
            else:
            	print("{:^80}".format("> Merging"))
                data = data.merge(df_tmp, how='outer', on=None)

            del processed_regions, df_tmp, df_tmp2
    opened_bw.close()

data.shape
data.columns
# Check number of lines in bed matches numer of rows of data
# !wc -l *siNipbl_log2Ratio*

# Missing Data as zero (np.ma.masked_array)
np.sum(data.iloc[:, 2:].values == 0)

# write to table
data.to_csv("data_log2scores_siNIPBL.tsv", sep="\t", header=True)

# Select columns with socres and .bed regions filename
colIx = [True if re.findall(r".+[.]bw_score", i)
         else False for i in data.columns.values]

cols = data.columns[colIx]
locIx = np.append(data.columns[colIx], "regions")

toPlot = pd.melt(data.loc[:, locIx], id_vars="regions",
                 value_vars=cols.values,
                 var_name="Mutant", value_name="Score")

toPlot.Mutant = [re.sub(r"_log2ratio.+", "", i) for i in toPlot.Mutant]
toPlot.regions = [re.sub(r"_log2Ratio.+", "", i) for i in toPlot.regions]

# numer obs and median
dfStats = toPlot.groupby(['Mutant', 'regions']).agg(['count', 'median'])

fig, ax = plt.subplots()
ax.axhline(y=0, color="red", ls='dashed')

ax = sns.boxplot(
    data=toPlot,
    x="Mutant",
    y="Score",
    hue="regions",
    notch=True)

ax = sns.stripplot(data=toPlot, x="Mutant", y="Score", hue="regions",
                   size=0.5, alpha=0.8, dodge=True)
# Add it to the plot
# pos = range(dfStats.shape[0])
# for tick,label in zip(pos,ax.get_xticklabels()):
# 	ax.text(pos[tick],
#             dfStats["Score","median"].values[tick] + 0.03,
#             dfStats["Score","count"].values[tick],
#             horizontalalignment='center',
#             size='x-small',
#             color='w',
#             weight='semibold')

# Remove repeated legends
hand, labl = ax.get_legend_handles_labels()
handout = []
lablout = []
for h, l in zip(hand, labl):
    if l not in lablout:
        lablout.append(l)
        handout.append(h)
ax.legend(handout, lablout)
ax.legend().remove()
fig.legend(handout, lablout)

# Format labels
ax.set_xticklabels(ax.get_xticklabels(), rotation=30)
ax.set_ylabel("Log2 Fold Change")
ax.set_xlabel("")

# save figure
plt.savefig(wd + '/score_siNIPBL@positions_ctcf-no_ctcf.png', res=300)
