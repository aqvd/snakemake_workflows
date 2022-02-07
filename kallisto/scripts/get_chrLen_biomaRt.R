# Alternative to get chrlens in R

library(GenomicFeatures)
getChromInfoFromBiomart(dataset = 'mmusculus_gene_ensembl',
                        host = "http://may2012.archive.ensembl.org") %>%
    dplyr::arrange(chrom) %>%
    head()