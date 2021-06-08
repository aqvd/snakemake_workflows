###########################################
#####    Find SA1/2 Only and Common   #####
###########################################

wd="/Users/aqo/Desktop/MCF10A/siNipbl/Peaks/merged_replicates_callpeak"

res_dir="${wd}/SAonly-common_allPeaks"
cd ${res_dir}

SA1_cont="${wd}/SA1_siC_merged_peaks.narrowPeak"
SA1_siN="${wd}/SA1_siNipbl_merged_peaks.narrowPeak"

SA2_cont="${wd}/SA2_siC_merged_peaks.narrowPeak"
SA2_siN="${wd}/SA2_siNipbl_merged_peaks.narrowPeak"

SMC1_cont="${wd}/SMC1_siC_merged_peaks.narrowPeak"
SMC1_siN="${wd}/SMC1_siNipbl_merged_peaks.narrowPeak"

wc -l ${wd}/{SA1,SA2,SMC1}*peaks.narrowPeak


## mergePraks compare peak files. Also split 1st line of mergePeaks output
## into command line and header line
mergePeaks -venn ${res_dir}/SA1_SA2_SMC1_siC_mergePeaks_venn.txt -d given ${SA1_cont} ${SA2_cont} ${SMC1_cont} | \
	gsed -E '1s/\)\tchr/\)\nchr/' 1> ${res_dir}/Cohesin_siC_mergePeaks.txt

mergePeaks -venn ${res_dir}/SA1_SA2_SMC1_siNipbl_mergePeaks_venn.txt -d given ${SA1_siN} ${SA2_siN} ${SMC1_siN} | \
	gsed -E '1s/\)\tchr/\)\nchr/' 1> ${res_dir}/Cohesin_siNipbl_mergePeaks.txt

## === Go to R to process mergePeaks out. venn_diagrams_cohesinOverlaps.R== #
## Concatenate to include STAG peaks that do not overlap with SMC1
cat SA1_siC.bed SA1_siC@SMC1_siC.bed > SA1_only_siC_all.bed
cat SA2_siC.bed SA2_siC@SMC1_siC.bed > SA2_only_siC_all.bed
cat SA1_siC@SA2_siC@SMC1_siC.bed SA1_siC@SA2_siC.bed > SA_common_siC_all.bed

cat SA1_siNipbl.bed SA1_siNipbl@SMC1_siNipbl.bed > SA1_only_siNipbl_all.bed
cat SA2_siNipbl.bed SA2_siNipbl@SMC1_siNipbl.bed > SA2_only_siNipbl_all.bed
cat SA1_siNipbl@SA2_siNipbl@SMC1_siNipbl.bed SA1_siNipbl@SA2_siNipbl.bed > SA_common_siNipbl_all.bed

ls *all.bed

## ========= Copy to cluster and plot heatmaps ========
## SA2_only_siNipbl are 1207 interesting positions. Use k means with respect to
## 27ac and 27me3 to separate into 8 clusters.

## Cluster 1 is unspeciffic, cluster7 could be those chrY regions. Check.  
cluster_bed="heatmap-CTCF_siNip_siSA_ring1_H3K27acMe3_zmym2@SA2only-siNipbl_NOTintersectedSMC1_Kmeans.bed"
head ${cluster_bed}

awk -F $'\t' 'BEGIN{OFS="\t"}; $13=="cluster_7" {print $1,$2,$3}' ${cluster_bed}
## Cluster 7 are those rubbish regions in chrY. Get rid of it too.

## Clusters 2 to 5 are K27ac high. Cluster 6 has a bit lower signal.
awk -F $'\t' 'BEGIN{OFS="\t"}; $13~/cluster_[2-5]/ {print $1,$2,$3,$13}' ${cluster_bed} > ${wd}/SA2_only_siNipbl_H3K27high.bed

awk -F $'\t' 'BEGIN{OFS="\t"}; $13~/cluster_6/ {print $1,$2,$3,$13}' ${cluster_bed} > ${wd}/SA2_only_siNipbl_H3K27mid.bed

