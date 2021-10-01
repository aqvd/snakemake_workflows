from snakemake.io import glob_wildcards, expand
import os
import pandas as pd
import numpy as np
import glob
import re
import functools

def field_from_sample(Sample, field):
    # ?P<group_name> for clarity
    match = re.match(
        r"(?P<readID>\w+)_(?P<readLetter>\w+)_(?P<readRun>\w+)_(?P<readLane>\w+)", Sample)

    return match.group(field)

GENOME_IX_PREFIX_DICT = {
    "mm9": "",
    "mm10": "",
    "hg19": "",
    "grch37": "",
    "hg38": "/home/aquevedo/resources/references/GRCh38/Verily_decoy/bwa-0.7.12/GRCh38_Verily_v1.genome.fa",
    "grch38": "/home/aquevedo/resources/references/GRCh38/Verily_decoy/bwa-0.7.12/GRCh38_Verily_v1.genome.fa",
    "-": ""}

REF_FASTA_DICT = {
    "hg38": "/home/aquevedo/resources/references/GRCh38/Verily_decoy/GRCh38_Verily_v1.genome.fa",
    "grch38": "/home/aquevedo/resources/references/GRCh38/Verily_decoy/GRCh38_Verily_v1.genome.fa",
    "mm10": ""
}

DBSNP_DICT = {
    "hg38": "/data_genome1/SharedSoftware/GATK/resources_hg38/dbsnp_146.hg38.vcf",
    "grch38": "/data_genome1/SharedSoftware/GATK/resources_hg38/dbsnp_146.hg38.vcf",
    "mm10": "-"
}

GNOMAD_SITES_DICT = {
    "hg38": "/data_genome1/SharedSoftware/GATK/resources_Mutect2/af-only-gnomad.hg38.vcf.gz",
    "grch38": "/data_genome1/SharedSoftware/GATK/resources_Mutect2/af-only-gnomad.hg38.vcf.gz",
    "mm10": "-"
}

GOLD_INDELS_DICT = {
    "hg38": "/data_genome1/SharedSoftware/GATK/resources_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
    "grch38": "/data_genome1/SharedSoftware/GATK/resources_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
    "mm10": "-"
}

PON_DICT = {
    "hg38": "/home/aquevedo/resources/GATK_public_data/1000g_pon.hg38.vcf.gz.tbi",
    "grch38": "/home/aquevedo/resources/GATK_public_data/1000g_pon.hg38.vcf.gz.tbi",
    "mm10": "-"
}

REGIONS_MUTECT_DIC = {
    "hg38": "/data_genome1/References/AgilentSureSelect/Human_Exome_V6_UTR/hg38/S07604624_Covered_fixed.bed",
    "grch38": "/data_genome1/References/AgilentSureSelect/Human_Exome_V6_UTR/hg38/S07604624_Covered_fixed.bed",
    "mm10": "-"
}

# Number chromosomes not including X (Y is discarded)
NCHR_DIC = {
    "hg38": 22,
    "hg19": 22,
    "grch38": 22,
    "grch37": 22,
    "mm10": 19,
    "mm9": 19
}

data = pd.read_csv(
    "/Users/aqo/Desktop/snakemake_workflows/test_WGS/table_test_GATK.txt", sep="\t")

data.columns
# Remove .fastq.gz to use basename with expand() in rule "all"
data["R1Basename"] = [f.replace(".fastq.gz", "") for f in data["R1"]]
data["R2Basename"] = [f.replace(".fastq.gz", "") for f in data["R2"]]

# Alignment
data["IxPrefPath"] = [GENOME_IX_PREFIX_DICT[i] for i in data.Genome]

data["ReadID"] = [field_from_sample(sample, "readID")
                  for sample in data.Samples]
data["ReadLetter"] = [field_from_sample(
    sample, "readLetter") for sample in data.Samples]
data["ReadRun"] = [field_from_sample(sample, "readRun")
                   for sample in data.Samples]
data["ReadLane"] = [field_from_sample(
    sample, "readLane") for sample in data.Samples]

# Base quality score recalibration parame
data["RefFASTA"] = [REF_FASTA_DICT[i] for i in data.Genome]
data["GoldIndels"] = [GOLD_INDELS_DICT[i] for i in data.Genome]
data["DbSNP"] = [DBSNP_DICT[i] for i in data.Genome]

# Mutec2 resources
data["PanelNormals"] = [PON_DICT[i] for i in data.Genome]
data["NumChr"] = [NCHR_DIC[i] for i in data.Genome]
data["RegionsMutect"] = [REGIONS_MUTECT_DIC[i] for i in data.Genome]


# TESTS 
def get_column_df(df , column, filt_out, **kwags):
    """
    Get column, eg TumorSample or NormalSample, for individual.
    Optional filter value for the results
    """
    # Filter using **kwags
    row_filters = []
    for k,v in kwags.items():
        row_filters.append('({0} == "{1}")'.format(k, v))

    query_expr = " & ".join(row_filters)

    # query expression and Select column with .loc[:, column]
    res = df.query(query_expr).loc[:, column]

    # Filter values of that column
    res = np.array(list(filter(lambda x: x != filt_out, res)))
    return(res)

    
a = get_column_df(data, "Samples", "", Individual = "indivA", IsControl = "no")
a

def repeat_argument(argument, value_list, **kwags):
    """
    sometimes a argument (eg: -I <file>) must be specified multiple times,
    this function creates the appropriate concatenation.
    Include whitespace in argument (eg "-I ") if needed
    """

    argument_l = [str(argument) + str(v) + " " for v in value_list]
    res = "".join(argument_l)

    return(res)

print(repeat_argument("-I ", 
      get_column_df(data, "Samples", "", Individual = "indivA", IsControl = "no")))

print(repeat_argument("-I ", 
      get_column_df(data, "Samples", "", Individual = "indivB", IsControl = "no")))

def expand_argument(path, string_to_expand,
                    df, column, filt_out,
                    argument, **kwags):
    """
    generates repeated argument value pairs (eg -I path/file_1.txt -I path/file_2.txt)
    using expand function from snakemake to create the paths.

    string_to_expand must be '{expansion}_bla_bla_bla.xxx'
    """
    if not path.endswith("/"):
        path = path + "/"

    # get the values from df
    to_expand = get_column_df(df, column, filt_out, **kwags)

    # expand them to file paths
    arg_value_l = expand(path + string_to_expand, expansion= to_expand)

    # Join using the argument
    res = repeat_argument(argument, arg_value_l)

    return(res)

DATADIR="/home/aquevedo/projects/test_GATKsnake/data/"
RESDIR="/home/aquevedo/projects/test_GATKsnake/res/"

get_column_df(data, "Samples", "", Individual = data.Individual[1], IsControl = "yes")[0]

expand_argument(DATADIR + "align", "{expansion}_rg_dedup_recal.bam",
     data, "Samples", "", "-I ", Individual = "indivA", IsControl = "no")

expand(DATADIR + "align/{sample}_rg_dedup_recal.bam",
                   sample = get_column_df(data, "Samples", "", 
                                          Individual = data.Individual.values[0],
                                          IsControl = "no"))
expand(RESDIR + "variants/{indiv}_somatic.vcf.gz", indiv = data.Individual.unique())


