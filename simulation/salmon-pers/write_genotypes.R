library(Biostrings)
library(tidyverse)

samples <- readLines("../samples.txt")

walk(samples, ~dir.create(file.path("quant", .)))

pers_transcripts <- "../../indices/personalize_transcripts/personalized_transcripts.fa" %>%
    readDNAStringSet()

individual_indices <- 
    "../../indices/personalize_transcripts/hla_genotypes_index.tsv" %>%
    read_tsv() %>%
    filter(sampleid %in% samples) %>%
    distinct(sampleid, gene_name, i, .keep_all = TRUE) %>%
    distinct(sampleid, allele) %>%
    split(.$sampleid) %>%
    map("allele") %>%
    map(unique) %>%
    map(~map(., ~grep(., names(pers_transcripts), value=TRUE, fixed = TRUE))) %>%
    map(unlist) %>%
    map(~pers_transcripts[.])

walk(samples, ~writeXStringSet(individual_indices[[.]], sprintf("quant/%s/hla.fa", .)))
