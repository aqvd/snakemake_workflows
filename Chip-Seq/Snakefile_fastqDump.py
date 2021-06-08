import os
import pandas as pd
import glob

## Directories must end with "/" !!!!!
FASTQDIR="/storage/scratch01/users/aquevedo/fastq/A673_chip_reanalysis/"
ACC_TABLE="/home/aquevedo/SRA_RunTables/SraRunTable_A673Chip_reanalysis.csv"

## Function to create directories unless they already exist
def tryMkdir(path):
	try:
		os.makedirs(path,exist_ok=False) # raise error if exists
	except FileExistsError:
		pass

tryMkdir(FASTQDIR)

data = pd.read_csv(ACC_TABLE,sep=",")

rule faster_download:
	input:
		expand(FASTQDIR + "{SRR}_faster.finished", SRR=data.Run.unique())

rule safe_download:
	input:
		expand(FASTQDIR + "{SRR}_prefetch.finished", SRR=data.Run.unique())

rule fastERq_dump:
	output:
		FASTQDIR + "{SRR}_faster.finished", # created only if all succeded on time
	params:
		srr=lambda wildcards: wildcards.SRR,
		Dir=FASTQDIR,
		prot=lambda wildcards: data.Protein[data.Run == wildcards.SRR].values[0],
		cond=lambda wildcards: data.Condition[data.Run == wildcards.SRR].values[0],
		rep=lambda wildcards: data.Rep[data.Run == wildcards.SRR].values[0]
	log:
		FASTQDIR + 'log/fastqDump/{SRR}.log'
	threads: 5
	shell:
		'fasterq-dump --log-level 5 --threads {threads} --temp {params.Dir} \
		--outfile {params.prot}_{params.cond}_{params.rep}_{params.srr}.fastq \
		--outdir {params.Dir} {params.srr} |& tee {log} && \
		gzip {params.Dir}{params.prot}_{params.cond}_{params.rep}_{params.srr}.fastq |& tee -a {log} && \
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
	threads: 5
	shell:
		'prefetch {params.srr} |& tee {log} && \
		fastq-dump --log-level 6 --outdir {params.Dir} {params.srr} |& tee -a {log} && \
		mv {params.Dir}{params.srr}.fastq \
		{params.Dir}{params.prot}_{params.cond}_{params.rep}_{params.srr}.fastq |& tee -a {log} && \
		gzip {params.Dir}{params.prot}_{params.cond}_{params.rep}_{params.srr}.fastq |& tee -a {log} && \
		touch {output} |& tee -a {log}'

	


