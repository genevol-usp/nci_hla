devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

samples <- readLines("../data/sample_ids.txt")

quants_noWin <-
    file.path("./quantifications_noWinner", samples, "quant.sf") %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "subject") %>%
    mutate(locus = imgt_to_gname(Name),
	   allele = sub("IMGT_", "", Name)) %>%
    select(subject, locus, allele, est_counts = NumReads, tpm = TPM)

typings_2nd <- quants_noWin %>%
    group_by(subject, locus) %>%
    slice(which.max(est_counts)) %>%
    ungroup() %>%
    filter(est_counts > 0)

typings_1st <- 
    read_tsv("./quantifications_topAlleles/processed_imgt_quants.tsv") %>%
    group_by(subject, locus) %>%
    slice(which.max(est_counts)) %>%
    ungroup()

typings_df <- bind_rows(typings_1st, typings_2nd) %>%
    arrange(subject, locus) %>%
    hla_genotype_dt(th = 0.01) %>%
    as_tibble()

calls <- typings_df %>%
    select(subject, locus, allele)

goldstd <- read_tsv("../../phase1/data/nci_expression.tsv") %>%
    select(subject, locus, allele)

calc_genotyping_accuracy(calls, goldstd) %>%
    write_tsv("./genotyping_concordance.tsv")

typings_df %>% 
    select(subject, locus, allele) %>%
    mutate(allele = paste0("IMGT_", allele)) %>%
    write_tsv("./genotype_calls.tsv")
