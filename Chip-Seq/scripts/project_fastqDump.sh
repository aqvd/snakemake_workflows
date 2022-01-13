#!/bin/bash

Script to create per project Snakefile-fastqDump

input:  
		projectID
		projectS directory
		fastq directory
		config table path
		snakemake-fastqDump path

output: 
		A modified snakemake-fastqDump in project directory
			- filename includes projectID
			- 