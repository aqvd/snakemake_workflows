snakemake -j1 -p -n --dag -s ~/snakemake_workflows/WGS/Snakefile_GATK.py \
--rerun-incomplete --forceall all


fil="/home/aquevedo/projects/test_WGS/res/variants/mouse2_somatic_filtered.vcf.gz"
unf="/home/aquevedo/projects/test_WGS/res/variants/mouse2_somatic_unfilt.vcf.gz"

cd $(dirname $fil)

bcftools view -H $fil |grep "PASS" | wc -l
bcftools view -H $unf |grep "PASS" | wc -l

bcftools view -H $fil | head

bcftools view -H $fil | cut -f7 | tr ";" "\n" |sort | uniq -c