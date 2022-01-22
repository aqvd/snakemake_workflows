#!/bin/bash

# Alvaro Quevedo
# Jan 2022
#
# >>    Configure the /home/aquevedo/resources directory   <<
#
# Copy, link, or download all the resources needed for variant calling
# both for human and mouse, according to GATK 4.2.0 best pracitces.
#
# Some data is already present in Genome1, like sequences and genome indexes.
# Other must be adquired from google cloud. interesting links:
#   - https://console.cloud.google.com/storage/browser/gatk-best-practices
#   - https://console.cloud.google.com/storage/browser/genomics-public-data
#

# Executables used
GATK="/data_genome1/SharedSoftware/GATK/gatk-4.1.9.0/gatk"
GSUTIL="/home/aquevedo/software/gsutil/gsutil"

# Create resource boundles for mouse and human genomes in server
basedir="/data_genome2/aquevedo/resources"
mkdir -p ${basedir}
cd ${basedir}

# ======= 1 make directories =====
parallel -k mkdir -p ${basedir}/genomes/{1}/{2}/{3} ::: "human" ::: "GRCh37-hg19" "GRCh38-hg38" ::: "sequences" "indexes" "variation" "annotation"
parallel -k mkdir -p ${basedir}/genomes/{1}/{2}/{3} ::: "mouse" ::: "GRCm37-mm9" "GRCm38-mm10" ::: "sequences" "indexes" "variation" "annotation"


## ============================================ ##
# ================= HUMAN hg38 =================
## ============================================ ##


hg38_dir="${basedir}/genomes/human/GRCh38-hg38"
cd ${hg38_dir}

# >> 2.1 Link sequences
# -- Already procssed sequences in genome 1
genome_dir="/data_genome1/References/Human/Sequences/Genome"
assembly="GRCh38_no_alt_plus_hs38d1_Verily"
ll ${genome_dir}/${assembly}/*

target_dir="${hg38_dir}/sequences/${assembly}"
mkdir -p ${target_dir}

ln -v -s --target-directory=${target_dir} ${genome_dir}/${assembly}/*

ll ${target_dir}

## -- Sequences processed by myself
target_dir="${hg38_dir}/sequences/GRCh38.p12-standardAndMito"
mkdir -p ${target_dir}

mv /data/aquevedo/genomes/human/GRCh38* $target_dir
ll -h ${target_dir}

ll ${hg38_dir}/indexes
mv ${target_dir}/GRCh38/*_bwa.* $hg38_dir/indexes
mv ${target_dir}/GRCh38/*.bt2 $hg38_dir/indexes
mv ${target_dir}/GRCh38/*.ht2 $hg38_dir/indexes
mv ${target_dir}/GRCh38/*build* $hg38_dir/indexes

mv ${target_dir}/GRCh38/*.gtf $hg38_dir/annotation
mv ${target_dir}/GRCh38/* ${target_dir}

mkdir -p $hg38_dir/indexes/GRCh38.p12-standardAndMito
mv $hg38_dir/indexes/* $hg38_dir/indexes/GRCh38.p12-standardAndMito

# >> 2.2 Link indexes present in "/data_genome1/References/GenomeIndices/Human/GRCh38_no_alt_plus_hs38d1" 
genome_dir="/data_genome1/References/GenomeIndices/Human"
assembly="GRCh38_no_alt_plus_hs38d1"
ll ${genome_dir}/${assembly}/*

target_dir="${hg38_dir}/indexes/${assembly}"
mkdir -p ${target_dir}

ln -v -s --target-directory=${target_dir} ${genome_dir}/${assembly}/*

ll ${target_dir}

# ======== 3. Extra resources needed for variant calling as GATK best practices ======
#	3.1- GNOmad
#	3.2- Panel Of Normals

# Use gsutil to acces data inn google cloud

target_dir="${hg38_dir}/variation"
cd $target_dir

$GSUTIL ls gs://gatk-best-practices/somatic-hg38/ | grep "vcf.gz" | parallel -j 6 ${GSUTIL} copy {} ${target_dir}/
ll -h ${target_dir}

bcftools view -H --threads 10 /data/aquevedo/resources/genomes/human/GRCh38-hg38/variation/af-only-gnomad.hg38.vcf.gz \
|head -n 100000000 | cut -f1 | sort | uniq -c


## ============================================ ##
# ================= MOUSE hg38 =================
## ============================================ ##


# Mouse lack GNOMAD resource, but SNPdb is present and could be used
# Hilmar did something similar as presented in this post and files are here
db_dir="/data_genome1/References/MusMusculus/Variation/dbSNP"

# check here https://www.biostars.org/p/338288/ and top rated answer linking to
# https://github.com/igordot/genomics/blob/master/workflows/gatk-mouse-mm10.md

mm10_dir="${basedir}/genomes/mouse/GRCm38-mm10"
cd $mm10_dir

# Download files parallely
# get list of files with curl -l and save to file, pipe is not working as expected
curl -l ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/ > download.files.tmp

awk '{print "ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/"$0}' download.files.tmp \
|grep -E "chr_([0-9XYM]{1,2}T?)\.+vcf.gz$" \
|parallel -j 10 wget --no-parent --no-directories --directory-prefix="${mm10_dir}/variation" {}

# Check version of dbSNP
zcat ${mm10_dir}/variation/vcf_chr_10.vcf.gz | head -n 20 #version 150
mkdir -p ${mm10_dir}/variation/dbSNP_v150

cd ${mm10_dir}/variation

# Add chr to CHROM names
for vcf in $(ls -1 vcf_chr_*.vcf.gz) ; do
  vcf_new=${vcf/.vcf.gz/.vcf}
  echo $vcf
  zcat $vcf | sed 's/^\([0-9XYMT]+\)/chr\1/' > $vcf_new
  rm -fv $vcf
done

ls -1 *.vcf | parallel bgzip {}
ls -1 *.vcf.gz | parallel bcftools index {}

bcftools merge *.vcf.gz | bcftools view --threads 6 -Oz > dbSNP_v150_all-chr-UCSC-GRCm38.vcf.gz

rm vcf*chr*vcf*
mv dbSNP_v150_all-chr-UCSC-GRCm38* dbSNP_v150/

# Add "chr" as it was not added with the previuos loop `
zcat dbSNP_v150/dbSNP_v150_all-chr-UCSC-GRCm38.vcf.gz \
| awk '{
	if ($0 !~ /^#/){
		print "chr"$0;
	} else {
		r = gensub( /(\#\#contig=<ID=)(.+)\>/, "\\1chr\\2", "g", $0);
		print r;}
}' \
|bcftools sort -m 4G -T ${mm10_dir} -Oz -o dbSNP_v150_all-chr-UCSCfix-GRCm38.vcf.gz

mv dbSNP_v150_all-chr-UCSCfix-GRCm38.vcf.gz ${mm10_dir}/variation/dbSNP_v150

zcat ${mm10_dir}/variation/dbSNP_v150/dbSNP_v150_all-chr-UCSCfix-GRCm38.vcf.gz \
|grep -E "^#" > ${mm10_dir}/variation/dbSNP_v150/dbSNP_v150_all-chr-UCSCfix-GRCm38.sorted &&
zcat ${mm10_dir}/variation/dbSNP_v150/dbSNP_v150_all-chr-UCSCfix-GRCm38.vcf.gz | grep -Ev "^#" \
| sort -V -k1,1 -k2,2n >> ${mm10_dir}/variation/dbSNP_v150/dbSNP_v150_all-chr-UCSCfix-GRCm38.sorted && \
| bgzip ${mm10_dir}/variation/dbSNP_v150/dbSNP_v150_all-chr-UCSCfix-GRCm38.sorted

mv ${mm10_dir}/variation/dbSNP_v150/dbSNP_v150_all-chr-UCSCfix-GRCm38.sorted.gz ${mm10_dir}/variation/dbSNP_v150/dbSNP_v150_all-chr-UCSC-GRCm38.vcf.gz

$GATK IndexFeatureFile -I ${mm10_dir}/variation/dbSNP_v150/dbSNP_v150_all-chr-UCSC-GRCm38.vcf.gz

# Check if AF can be added to VCF file like posted here:
# https://gatk.broadinstitute.org/hc/en-us/community/posts/360078092732-Common-germline-variants-sites-VCF-for-GetpileupSummaries-for-mouse-data

wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/
less README

# Download C57BL-6NJ to see the INFO fields present in those vcfs
curl -s -l ftp://ftp-mouse.sanger.ac.uk/current_snps/
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz


bcftools view -h C57BL_6NJ.mgp.v5.indels.dbSNP142.normed.vcf.gz | grep "INFO"
bcftools view -H C57BL_6NJ.mgp.v5.indels.dbSNP142.normed.vcf.gz \
|head -n 3 \
| 

bcftools +fill-tags C57BL_6NJ.mgp.v5.indels.dbSNP142.normed.vcf.gz -- -t AF | \
grep -Eo ";AF=[0-9.]+" | 
grep -Eo "[0-9.]+" |
sort | uniq -c

# solution here: https://www.biostars.org/p/180894/
# find bcftools
which bcftools | basename
# set environmental variable to use plugins
export BCFTOOLS_PLUGINS=/data_genome1/SharedSoftware/samtools_htslib/download/bcftools-1.12/plugins

# test plugins
bcftools view -H C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf.gz | head -5
zcat C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf.gz | head -n 100 | \

bcftools +fill-tags C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf.gz -- -t AF | \
grep -v "^#" | \
cut -f 8 | \
grep -Eo "AF=([0-9])+" |
sed -E "s/AF=([0-9])+/\1/" | \
sort | uniq -c



# Install latest bcftools and samtools to be able to use plugins
git clone git://github.com/samtools/htslib.git
git clone git://github.com/samtools/bcftools.git
cd bcftools; ./configure --prefix=/home/aquevedo/software/ make
sudo make install










tabix \
  ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz \
  2:3000000-3020000

# generate parameter string containing all VCF files

# create new vcf with hromosomes in UCSC style, "chr[0-9XYMT]"
#	- v146
# cat <(zcat ${mm10_dir}/variation/dbSNP_v146/dbSNP_all.chr.vcf.gz | head -n 30 | egrep "^#") \
# <(zcat ${mm10_dir}/variation/dbSNP_v146/dbSNP_all.chr.vcf.gz | egrep -v "^#" |awk '{print "chr"$0;}') \
# | bgzip > ${mm10_dir}/variation/dbSNP_v146/dbSNP_chr-all-UCSC.vcf.gz

# zcat ${mm10_dir}/variation/dbSNP_v146/dbSNP_chr-all-UCSC.vcf.gz | head -n 30 
# $GATK IndexFeatureFile -I ${mm10_dir}/variation/dbSNP_v146/dbSNP_chr-all-UCSC.vcf.gz 


## >> INDELS, taken form mouse genomes project.
# To browse the mgp ftp server use "curl -l ftp://ftp-mouse.sanger.ac.uk"

target_dir="${mm10_dir}/variation/indels"
mkdir -p ${mm10_dir}/variation/indels

cd ${target_dir}

wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz \
-O ${target_dir}/mgp.v5.indels.vcf.gz

# adjust header
zcat ${target_dir}/mgp.v5.indels.vcf.gz | head -1000 | grep "^#" | cut -f 1-8 \
| grep -v "#contig" | grep -v "#source" \
> ${target_dir}/mgp.v5.indels.pass.chr.vcf
# keep only passing and adjust chromosome name
zcat ${target_dir}/mgp.v5.indels.vcf.gz | grep -v "^#" | cut -f 1-8 \
| grep -w "PASS" | sed -E 's/^([0-9XYMT]+)/chr\1/' | cut -f 1 | uniq -c 
>> ${target_dir}/mgp.v5.indels.pass.chr.vcf

# compress the vcf and index
bcftools view -Oz --threads 6 ${target_dir}/mgp.v5.indels.pass.chr.vcf > ${target_dir}/mgp.v5.indels.pass.chr.vcf.gz 

$GATK IndexFeatureFile -I ${target_dir}/mgp.v5.indels.pass.chr.vcf.gz
zcat ${target_dir}/mgp.v5.indels.pass.chr.vcf.gz | head -n 50

# remove uncompressed and no filtered .vcf
rm ${target_dir}/mgp.v5.indels.pass.chr.vcf ${target_dir}/mgp.v5.indels.vcf.gz




