# Generate file equivalent to GNOMAD resource for mm10

# gnomad-af resource for GATK requires Allele Frequency field (AF)
# this is computed as AF=AC/AN 
# (https://gatk.broadinstitute.org/hc/en-us/community/posts/360069566651-Mutect2-germline-resources-for-mm10)
#
# So, from mouse genomes project, one can extract the .vcf for diffrent mouse strains (healthy)
# and estimate AF (if not already present) using that formula
#
dir="/data/aquevedo/resources/genomes/mouse/GRCm38-mm10/variation/mgp_snps_indels"
mkdir -p $dir
cd $dir

# list ftp directory and go to the latest release (Apr20)
curl -l ftp://ftp-mouse.sanger.ac.uk/REL-2004-v7-SNPs_Indels/

# Use tabix to query a region and see if AF, or at least AC and AN fields, are present
tabix ftp://ftp-mouse.sanger.ac.uk/REL-2004-v7-SNPs_Indels/mgp_REL2005_snps_indels.vcf.gz 2:3000000-4000000 \
|head -n 3 \
|bcftools +fill-tags -- -t AF
|sed -En 's/.+(AC=[0-9]+;AN=[0-9]+);.+/\1/p' \

# In order not to download 39Gb of file, use curl to process while download
# Keep header, Column names line and fields 1-8 of .vcf
cat <(curl -s ftp://ftp-mouse.sanger.ac.uk/REL-2004-v7-SNPs_Indels/mgp_REL2005_snps_indels.vcf.gz \
	  | zcat \
	  | head -n 100 \
	  | grep -E "^##" ) \
	<(curl -s ftp://ftp-mouse.sanger.ac.uk/REL-2004-v7-SNPs_Indels/mgp_REL2005_snps_indels.vcf.gz \
		  | zcat \
		  | head -n 100 \
		  | grep -E "^#CHROM" \
		  | cut -f 1-8) \
    <(curl -s ftp://ftp-mouse.sanger.ac.uk/REL-2004-v7-SNPs_Indels/mgp_REL2005_snps_indels.vcf.gz \
    	| zcat \
    	| grep -Ev "^#" \
	    | grep "PASS" \
	    | cut -f 1-8 ) > mgp_REL2005_snps_indels_f1-8.vcf

gzip mgp_REL2005_snps_indels_f1-8.vcf
# I dont know why but bcftools plugin fill-tags is not working in this vcf,
# So compute manually using awk and then paste that column
zcat mgp_REL2005_snps_indels_f1-8.vcf.gz | head -n 80 

# Thrar is the ##INFO field for AF
bcftools +fill-tags mgp_REL2005_snps_indels_f1-8.vcf -- -t AF | head -n 80 | grep AF
# ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">

# Compute AF 
zcat mgp_REL2005_snps_indels_f1-8.vcf.gz \
|awk '\
{ if ($0 !~ /^#/ ) { 
	AC = gensub(/.+;AC=([0-9]+).+/, "\\1", "g", $8);
	AN = gensub(/.+;AN=([0-9]+).+/, "\\1", "g", $8);
	print "AF="AC/AN;
	}
}' \
> AF.txt

gunzip mgp_REL2005_snps_indels_f1-8.vcf.gz

# get header, add AF field
cat <(head  -n 80 mgp_REL2005_snps_indels_f1-8.vcf | grep "^##") \
 	<(echo -e '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">') \
 	<(head  -n 80 mgp_REL2005_snps_indels_f1-8.vcf | grep "^#CHROM") > mgp_REL2005_snps_indels_f1-8_AF.vcf

cat mgp_REL2005_snps_indels_f1-8_AF.vcf
# check same number of lines before pasting
grep -v "^#" mgp_REL2005_snps_indels_f1-8.vcf | wc -l 
wc -l AF.txt 


cat <(paste --delimiters="" \
 	  <(printf -- "AF=%0.s\n" {1..10}) \
	  <(head -n 10 AF.txt)) 

# add AF field
paste --delimiters=";" \
 		    <(grep -v "^#" mgp_REL2005_snps_indels_f1-8.vcf) \
 		    <(paste --delimiters="" \
			 	  <(printf -- "AF=%0.s\n" {1..92465298}) \
				  AF.txt \
			  ) >> mgp_REL2005_snps_indels_f1-8_AF.vcf

bcftools view --threads 10 -Oz mgp_REL2005_snps_indels_f1-8_AF.vcf > mgp_REL2005_snps_indels_f1-8_AF.vcf.gz

head -n 100 mgp_REL2005_snps_indels_f1-8_AF.vcf

# remove vcf without AF
# rm mgp_REL2005_snps_indels_f1-8.vcf mgp_REL2005_snps_indels_f1-8_AF.vcf

bcftools view mgp_REL2005_snps_indels_f1-8_AF.vcf.gz | head -n 100

# SNPS
cat <(bcftools view -h mgp_REL2005_snps_indels_f1-8_AF.vcf.gz) \
	<(bcftools view --threads 10 -Ov -H mgp_REL2005_snps_indels_f1-8_AF.vcf.gz \
		| grep -vE "INDEL;") \
|bcftools view --threads 10 -Oz - > mgp_REL2005_snps_f1-8_AF.vcf.gz &

bcftools view mgp_REL2005_snps_f1-8_AF.vcf.gz | head -n 100

# INDELS
cat <(bcftools view -h mgp_REL2005_snps_indels_f1-8_AF.vcf.gz) \
	<(bcftools view --threads 10 -Ov -H mgp_REL2005_snps_indels_f1-8_AF.vcf.gz \
		| grep -E "INDEL;") \
|bcftools view --threads 10 -Oz - > mgp_REL2005_indels_f1-8_AF.vcf.gz &

# Add chr to chr names
for f in mgp*_AF.vcf.gz; do
	l=$(bcftools view -H $f | wc -l | sed -Es 's/[[:space:]].+//g')
	printf -- "------ > %s < ------\n lines = %s \n\n" $f $l

	cat <(bcftools view -h $f) \
		<(paste --delimiters="" \
				<(printf -- "chr\n%0.s" $(seq 1 ${l}) ) \
				<(bcftools view --threads 10 -Ov -H $f ) \
		  ) \
	|bcftools view --threads 10 -Oz - > ${f/.vcf.gz/_chrUCSC.vcf.gz} 
done

# Sort VCFs by coordinate
ls -1 *_chrUCSC.vcf.gz \
|parallel -j 3 --tmpdir=/data/aquevedo/resources/ '\
	cat <(bcftools view -h -Ov {}) \
		<(bcftools view --threads 5 -H -Ov {} |sort -V -k1,1 -k2,2n) \
	|bcftools view --threads 5 -Oz - > {.}.sorted.gz '

rm *_chrUCSC.vcf.gz *_AF.vcf.gz

# index using tabix, also fix CONTIG header info after adding chr
ls -1 *_chrUCSC.vcf.sorted.gz | parallel --tmpdir=/data/aquevedo/resources/ tabix -p vcf {}
rm *_chrUCSC.vcf.sorted.gz.tbi

# I don't know why, but contig info appears two times. First time without "chr"
# Find which lines are thoser and delete
ls -1 *chrUCSC.vcf.sorted.gz \
|parallel --tmpdir=/data/aquevedo/resources/ \
'printf -- "\n------- > %s < --------\n" {} ; \
bcftools view -h {} | grep -nE "##contig=<ID=[0-9X]" | cut -d ":" -f 1 | tr "\n" ","'

# Delete lines 7 to 26
ls -1 *chrUCSC.vcf.sorted.gz \
|parallel --tmpdir=/data/aquevedo/resources/ \
'printf -- "\n------- > %s < --------\n" {} ; \
bcftools view -Ov --threads 10 {} \
| sed "7,26d" \
|bcftools view -Oz --threads 10 - > {.}.fix.gz' 

#check
ls -1 *chrUCSC.vcf.sorted.fix.gz \
|parallel --tmpdir=/data/aquevedo/resources/ \
'printf -- "\n------- > %s < --------\n" {} ; \
bcftools view -h {} | grep -nE "##contig=<ID="'

ls -1 *chrUCSC.vcf.sorted.fix.gz \
|parallel --tmpdir=/data/aquevedo/resources/ \
'mv {} "{= s/.fix.gz/.gz/ =}"'

# ====== Select Biallelic with >0.1 AF and Remove unnecessary fields ===== #

# - Script used for GATK team to build the resource
#     https://github.com/broadinstitute/gatk/blob/master/scripts/mutect2_wdl/mutect_resources.wdl
#
# - Discussion where 0.1 used as cutoff. Answer of developer of GATK
#	  https://gatk.broadinstitute.org/hc/en-us/community/posts/360072525211-Issues-with-Java-heap-space-in-CollectAllelicCounts
#
dir="/data/aquevedo/resources/genomes/mouse/GRCm38-mm10/variation/mgp_snps_indels"
mkdir -p $dir
cd $dir

GATK="/data_genome1/SharedSoftware/GATK/gatk-4.1.9.0/gatk"

ls -1 *chrUCSC.vcf.sorted.gz \
|parallel  --tmpdir=/data/aquevedo/resources/ \
$GATK --java-options \"-Xmx4096M -Djava.io.tmpdir=/data/aquevedo/resources/\" \
SelectVariants \
-V {} \
--restrict-alleles-to BIALLELIC \
-select \"AF \> 0.1\" \
-O "{= s/.vcf.sorted.gz/_biallelic.sorted.vcf.gz/ =}"

bcftools view -H mgp_REL2005_indels_f1-8_AF_chrUCSC_biallelic.sorted.vcf.gz | head -n 2

ls -1 *_biallelic.sorted.vcf.gz \
|parallel --tmpdir=/data/aquevedo/resources/ \
'cat <(bcftools view -h {}) \
	 <(paste -d $"\t" \
		    <(bcftools view --threads 6 -H {} | cut -f 1-7) \
		    <(bcftools view --threads 6 -H {} | cut -f 8 \
		    	| cut -d ";" -f 1,3) \
	  ) \
| bcftools view -Oz - > "{= s/_biallelic./_afOnly_bialellelic./ =}"' &

bcftools view -H  mgp_REL2005_snps_f1-8_AF_chrUCSC_afOnly_bialellelic.sorted.vcf.gz | head 
bcftools view -H mgp_REL2005_indels_f1-8_AF_chrUCSC.vcf.sorted.gz | cut -f1 | uniq -c
bcftools view -h  mgp_REL2005_snps_indels_f1-8_AF_chrUCSC_afOnly_bialellelic.sorted.vcf.gz

rm *_afOnly_bialellelic.*.tbi

# Index using GATK to be sure it wont crash in pipeline
ls -1 *_afOnly_bialellelic.* | parallel -j 3 --tmpdir="/data/aquevedo/resources/" \
/data_genome1/SharedSoftware/GATK/gatk-4.1.9.0/gatk \
--java-options \"-Xmx4096M -Djava.io.tmpdir=/data/aquevedo/resources/\" \
IndexFeatureFile -I {}

less 


