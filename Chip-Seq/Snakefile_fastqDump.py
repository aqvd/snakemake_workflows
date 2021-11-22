import os
import pandas as pd
import glob

## Directories must end with "/" !!!!!
FASTQDIR="/home/aquevedo/fastq/SraRunTable-GSE144833_SUIT2control_rna_exome/"
ACC_TABLE="/home/aquevedo/projects/SUIT2control_rna_exome/SraRunTable_SUIT2cont_GSE144833.csv"

##Function to create directories unless they already exist
def tryMkdir(path):
	try:
		os.makedirs(path,exist_ok=False) # raise error if exists
	except FileExistsError:
		pass

tryMkdir(FASTQDIR)

data = pd.read_csv(ACC_TABLE,sep=",")

rule faster_download:
	input:
		expand(FASTQDIR + "{SRR}_faster_gzip.finished", SRR=data.Run.unique())

rule safe_download:
	input:
		expand(FASTQDIR + "{SRR}_prefetch_gzip.finished", SRR=data.Run.unique())

rule fastERq_dump:
	output:
		FASTQDIR + "{SRR}_faster.finished", # created only if all succeded on time
	params:
		srr=lambda wildcards: wildcards.SRR,
		Dir=FASTQDIR
	log:
		FASTQDIR + 'log/fastqDump/{SRR}.log'
	threads: 4
	shell:
		'fasterq-dump --threads {threads} --temp {params.Dir} --split-3 --skip-technical \
		--outdir {params.Dir} {params.srr} |& tee {log} && \
		touch {output} |& tee -a {log}'

rule prefetch_fastq_dump:
	output:
		FASTQDIR + "{SRR}_prefetch.finished", # created only if all succeded on time
	params:
		srr=lambda wildcards: wildcards.SRR,
		Dir=FASTQDIR,
		prot=lambda wildcards: data.Protein[data.Run == wildcards.SRR].values[0],
		cond=lambda wildcards: data.Condition[data.Run == wildcards.SRR].values[0],
		rep=lambda wildcards: data.Rep[data.Run == wildcards.SRR].values[0]
	log:
		FASTQDIR + 'log/fastqDump/{SRR}.log'
	threads: 1
	shell:
		'prefetch {params.srr} |& tee {log} && \
		fastq-dump --log-level 6 --split-3 --skip-technical --dumpbase --clip \
		--outdir {params.Dir} {params.srr} |& tee -a {log} \
		touch {output} |& tee -a {log}'

rule pigz_and_rename_fq:
	input:
		FASTQDIR + "{SRR}_{method}.finished"
	output:
		FASTQDIR + "{SRR}_{method}_gzip.finished"
	params:
		srr=lambda wildcards: wildcards.SRR,
		Dir=FASTQDIR,
		prot=lambda wildcards: data.Protein[data.Run == wildcards.SRR].values[0],
		cond=lambda wildcards: data.Condition[data.Run == wildcards.SRR].values[0],
		rep=lambda wildcards: data.Rep[data.Run == wildcards.SRR].values[0]
	log:
		FASTQDIR + 'log/fastqDump/{SRR}_{method}_gzip.log'
	threads: 3
	shell:
		'scripts/name_split3.sh {params.Dir} {params.srr} {params.prot} \
								params.cond} {params.rep} {threads} |& tee -a {log}'




