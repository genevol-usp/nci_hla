library(tidyverse)

dir.create("genotypes")

fasta <- Biostrings::readDNAStringSet("./hladb/transcripts_MHC_HLAsupp.fa") %>%
    .[grepl("^IMGT", names(.))]

genos <- read_tsv("~/simulation_rnaseq/hlagenotypes.tsv") %>%
    mutate(locus = sub("HLA-", "", gene_name),
	   allele = paste0("IMGT_", allele)) %>%
    select(subject, locus, allele) %>%
    mutate(allele2 = ifelse(grepl("(:\\d+){3}", allele),
			    sub("(:\\d+$)", "", allele),
			    allele)) %>%
    mutate(final = case_when(allele %in% names(fasta) ~ allele,
			     !allele %in% names(fasta) & allele2 %in% names(fasta) ~ allele2,
			     TRUE ~ NA_character_)) %>%
    select(subject, locus, allele = final)

if (any(is.na(genos$allele))) stop("missing allele")

genos %>%
    split(.$subject) %>%
    walk(~write_tsv(., sprintf("./genotypes/%s.tsv", unique(.$subject))))

