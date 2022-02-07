params_index="/data/aquevedo/resources/genomes/human/GRCh38-hg38/indexes/ensembl_103/Homo_sapiens.GRCh38.cdna.all.release-103_k31.idx"
params_outdir="/data/aquevedo/projects/p53_ARID1A_DKO-GSE164179/kallisto/"
fqdir="/data/aquevedo/fastq/p53_ARID1A_DKO-GSE164179/"
gtf="/data/aquevedo/resources/genomes/human/GRCh38-hg38/annotation/Homo_sapiens.GRCh38.103.gtf.gz"
chr_len="/data/aquevedo/resources/genomes/human/GRCh38-hg38/sequences/ensembl_103/chrLens.txt"

kallisto quant \
--index="${params_index}" \
--output-dir="${params_outdir}" \
--bootstrap-samples=4 \
--threads=2 \
--genomebam \
--gtf "${gtf}" \
--chromosomes ${chr_len} \
<(zcat ${fqdir}TP53-ARID1A_DKO_R1_SRR13347775_R1.fastq.gz) \
<(zcat ${fqdir}TP53-ARID1A_DKO_R1_SRR13347775_R2.fastq.gz)



ll ${params_index}
ll ${params_outdir}

<(zcat ${fqdir}TP53-ARID1A_DKO_R1_SRR13347777_R1.fastq.gz) \
<(zcat ${fqdir}TP53-ARID1A_DKO_R1_SRR13347777_R2.fastq.gz) \
<(zcat ${fqdir}TP53-ARID1A_DKO_R1_SRR13347778_R1.fastq.gz) \
<(zcat ${fqdir}TP53-ARID1A_DKO_R1_SRR13347778_R2.fastq.gz) \
<(zcat ${fqdir}TP53-ARID1A_DKO_R1_SRR13347779_R1.fastq.gz) \
<(zcat ${fqdir}TP53-ARID1A_DKO_R1_SRR13347779_R2.fastq.gz) \
<(zcat ${fqdir}TP53-ARID1A_DKO_R1_SRR13347780_R1.fastq.gz) \
<(zcat ${fqdir}TP53-ARID1A_DKO_R1_SRR13347780_R2.fastq.gz) \
<(zcat ${fqdir}TP53-ARID1A_DKO_R1_SRR13347781_R1.fastq.gz) \
<(zcat ${fqdir}TP53-ARID1A_DKO_R1_SRR13347781_R2.fastq.gz) \
<(zcat ${fqdir}TP53-ARID1A_DKO_R1_SRR13347782_R1.fastq.gz) \
<(zcat ${fqdir}TP53-ARID1A_DKO_R1_SRR13347782_R2.fastq.gz) 



kallisto index -i transcripts.kidx transcripts.fasta.gz
kallisto quant         -i transcripts.kidx         -b 30         -o quant_out/         --genomebam         --gtf transcripts.gtf.gz         --chromosomes chrom.txt         reads_1.fastq.gz reads_2.fastq.gz


params_index="/data/aquevedo/resources/genomes/human/GRCh38-hg38/indexes/ensembl_103/Homo_sapiens.GRCh38.cdna.all.release-103_k31.idx"
params_outdir="/data/aquevedo/projects/p53_ARID1A_DKO-GSE164179/kallisto/"
fqdir="/data/aquevedo/fastq/p53_ARID1A_DKO-GSE164179/"
gtf="/data/aquevedo/resources/genomes/human/GRCh38-hg38/annotation/Homo_sapiens.GRCh38.103.gtf.gz"
chr_len="/data/aquevedo/resources/genomes/human/GRCh38-hg38/sequences/ensembl_103/chrLens.txt"

kallisto quant \
--index="Homo_sapiens.GRCh38.cdna.all.release-103_k31.idx" \
--output-dir="kallisto/" \
--bootstrap-samples=4 \
--threads=2 \
--genomebam \
--gtf "Homo_sapiens.GRCh38.103.gtf.gz" \
--chromosomes "chrLens.txt" \
<(zcat SRR13347775_R1.fastq.gz) \
<(zcat SRR13347775_R2.fastq.gz)


conda create -n kallisto
conda activate kallisto
mamba install -c conda-forge -c bioconda kallisto

kallisto quant -i transcripts.kidx \
-b 30 -o quant_out/ \
--genomebam --gtf transcripts.gtf.gz \
--chromosomes chrom.txt \
reads_1.fastq.gz reads_2.fastq.gz