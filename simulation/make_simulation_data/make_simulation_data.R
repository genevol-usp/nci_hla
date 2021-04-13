library(Biostrings)
library(hlaseqlib)
library(tidyverse)

# samples
nci_samples <- read_tsv("../../analysis/samples.txt", col_names = c("tp", "sampleid")) %>%
    filter(tp == 1)

# personalized transcripts
pers_transcripts <- "../../indices/personalize_transcripts/personalized_transcripts.fa" %>%
    readDNAStringSet()

# annotations
annots <- "/home/vitor/gencode/gencode.v37.primary_assembly.annotation.gtf" %>% 
    read_tsv(comment = "#", col_names = FALSE, col_types = "c-cii-c-c")

gene_annots <- annots %>%
    filter(X3 == "gene") %>%
    transmute(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
              gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
              strand = X7) %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C"))

# select protein-coding transcripts
hla_transcripts_annot <- annots %>%
    filter(X3 == "transcript") %>%
    transmute(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	      gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"), 
	      tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
	      tx_type = str_extract(X9, "(?<=transcript_type\\s\")[^\"]+")) %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C"))

hla_transcripts_protein <- hla_transcripts_annot %>%
    filter(tx_type == "protein_coding") %>%
    select(-tx_type)

## select transcripts which explain 90% of the expression in the real data
nci_expression <- "../../analysis/salmon/quant/%s_t1/quant.sf" %>%
    sprintf(nci_samples$sampleid) %>%
    setNames(nci_samples$sampleid) %>%
    map_df(read_tsv, .id = "sampleid") %>%
    group_by(sampleid) %>%
    mutate(scaled_counts = as.integer(round(NumReads/sum(NumReads) * 3e7))) %>%
    ungroup()

hla_expression <- inner_join(nci_expression, hla_transcripts_protein, by = c("Name" = "tx_id")) 

hla_selected_transcripts <- hla_expression %>%
    group_by(sampleid, gene_name) %>%
    mutate(avg_tpm = (TPM/sum(TPM)) * 100) %>%
    group_by(gene_name, tx_id = Name) %>%
    summarise(avg_tpm = mean(avg_tpm)) %>%
    ungroup() %>%
    arrange(gene_name, -avg_tpm) %>%
    group_by(gene_name) %>%
    mutate(cumexp = cumsum(avg_tpm),
	   i = cumsum(cumexp >= 90)) %>%
    ungroup() %>%
    filter(i <= 1L)

# HLA genotypes
hla_genotypes <- read_tsv("../../indices/personalize_transcripts/hla_genotypes_index.tsv") %>%
    filter(sampleid %in% hla_expression$sampleid) %>%
    distinct(sampleid, gene_name, i, .keep_all = TRUE)

## choose 50 individuals
set.seed(10)
simul_df <- hla_genotypes %>%
    group_by(sampleid) %>%
    nest() %>%
    ungroup() %>%
    arrange(sampleid) %>%
    sample_n(50) %>%
    unnest(c(data)) %>%
    arrange(sampleid)

simul_df %>%
    distinct(sampleid) %>%
    pull(sampleid) %>%
    writeLines("../samples.txt")

# hla index
hla_index <- 
    tibble(id = names(pers_transcripts)) %>%
    separate(id, c("tx_id", "allele"), sep = "_", remove = FALSE) %>%
    filter(allele %in% unique(simul_df$allele)) %>%
    pull(id) %>%
    pers_transcripts[.]

# make count matrix
genome_background_quants <- nci_expression %>%
    filter(sampleid == "66K00003") %>%
    select(tx_id = Name, readcount = scaled_counts) %>%
    filter(! tx_id %in% hla_transcripts_annot$tx_id)

hla_ground_truth <- hla_expression %>%
    filter(Name %in% unique(hla_selected_transcripts$tx_id)) %>%
    select(sampleid, gene_name, tx_id = Name, readcount = scaled_counts) %>%
    arrange(sampleid, gene_name, tx_id) %>%
    inner_join(simul_df, by = c("sampleid", "gene_name")) %>%
    group_by(sampleid, tx_id) %>%
    mutate(readcount = readcount/2L) %>%
    ungroup() %>%
    unite(tx_id, c("tx_id", "allele"), sep = "_") %>%
    group_by(sampleid, tx_id) %>%
    summarise(readcount = sum(readcount)) %>%
    ungroup()

final_quants <- hla_ground_truth %>%
    split(.$sampleid) %>%
    map(~select(., -sampleid)) %>%
    map_df(~bind_rows(genome_background_quants, .), .id = "sampleid") %>%
    group_by(sampleid) %>%
    mutate(readcount = as.integer(round(readcount/sum(readcount) * 3e7))) %>%
    ungroup()

pheno_matrix <- final_quants %>%
    pivot_wider(names_from = sampleid, values_from = readcount) %>%
    mutate_at(vars(-tx_id), ~replace_na(., 0))

# index for simulation
index_transcriptome <- "../../indices/gencode.transcripts.fa" %>%
    readDNAStringSet()

index <- c(index_transcriptome, hla_index)[pheno_matrix$tx_id]

# save output
write_tsv(select(pheno_matrix, -1), "phenotypes.tsv")
writeXStringSet(index, "simulation_index.fa")

# compute mean and sd of fragment length distribution
fld <- scan("../../analysis/salmon/quant/66K00003_t1/libParams/flenDist.txt")

fld_params <- tibble(mu = sum(1:length(fld) * fld)) %>%
    mutate(sigma = sqrt(sum( ( (1:length(fld) - mu)^2 ) * (fld) )))

write_tsv(fld_params, "fld_params.tsv")
