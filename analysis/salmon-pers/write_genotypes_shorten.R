library(Biostrings)
library(tidyverse)

samples <- read_tsv("../samples.txt", col_names = c("phase", "sampleid")) %>%
    mutate(id = sprintf("%s_t%s", sampleid, phase))

dir.create("quant-shorten")
walk(samples$id, ~dir.create(file.path("quant-shorten", .)))

pers_transcripts <- "../../indices/personalize_transcripts/personalized_transcripts_plusShort.fa" %>%
    readDNAStringSet()

individual_indices <- 
    "../../indices/personalize_transcripts/hla_genotypes_index.tsv" %>%
    read_tsv() %>%
    inner_join(samples, by = "sampleid") %>%
    distinct(id, gene_name, i, .keep_all = TRUE) %>%
    distinct(id, allele) %>%
    arrange(id, allele) %>%
    split(.$id) %>%
    map("allele") %>%
    map(unique) %>%
    map(~map(., ~grep(., names(pers_transcripts), value=TRUE, fixed = TRUE))) %>%
    map(unlist) %>%
    map(~pers_transcripts[.])

walk(samples$id, ~writeXStringSet(individual_indices[[.]], sprintf("quant-shorten/%s/hla.fa", .)))
