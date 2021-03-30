library(tidyverse)

hlagenos <- list.files("./hlapers/genotypes", full.names = TRUE) %>%
    map_df(read_tsv) %>%
    filter(locus %in% c("A", "B", "C")) %>%
    mutate(gene_name = paste0("HLA-", locus)) %>%
    select(gene_name, tx_id = allele) %>%
    distinct(gene_name, tx_id)

annots <- "~/gencode/gencode.v37.primary_assembly.annotation.gtf" %>%
    read_tsv(comment = "#", col_names = FALSE, col_types = "c-cii-c-c")

gene_annots <- annots %>%
    filter(X3 == "gene") %>%
    transmute(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	      gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"))

tx_annots <- annots %>%
    filter(X3 == "transcript") %>%
    transmute(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	      gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	      tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"))

hla_annots <- tx_annots %>%
    filter(gene_name %in% paste0("HLA-", c("A", "B", "C"))) %>% 
    bind_rows(hlagenos) %>%
    split(.$gene_name) %>%
    map_df(~fill(., gene_id, .direction = "updown"))
	    

samples <- read_tsv("./samples.txt", col_names = c("timepoint", "sampleid"))

sample_ids <- sprintf("%s_t%s", samples$sampleid, samples$timepoint)


###
salmon_t1 <- file.path("./salmon/quant", sample_ids, "quant.sf") %>%
    setNames(sample_ids) %>% .[grepl("_t1", .)] %>%
    map_df(~read_tsv(.) %>%
	   left_join(tx_annots, by = c("Name" = "tx_id")) %>%
	   group_by(gene_id, gene_name) %>%
	   summarise(tpm = sum(TPM)) %>%
	   ungroup(), .id = "sampleid")

write_rds(salmon_t1, "./plot_data/salmon_quants.rds")


###

salmon <- file.path("./salmon/quant", sample_ids, "quant.sf") %>%
    setNames(sample_ids) %>%
    map_df(read_tsv, .id = "sampleid") %>%
    group_by(sampleid) %>%
    mutate(cpm = (NumReads/sum(NumReads)) * 1e6) %>%
    ungroup() %>%
    inner_join(hla_annots, by = c("Name" = "tx_id")) %>%
    group_by(sampleid, gene_name) %>%
    summarise(cpm = sum(cpm)) %>%
    ungroup() %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_")

hlapers <- 
    sprintf("./hlapers/quant/%s_quant/quant.sf", sample_ids) %>%
    setNames(sample_ids) %>%
    map_df(read_tsv, .id = "sampleid") %>%
    group_by(sampleid) %>%
    mutate(cpm = (NumReads/sum(NumReads)) * 1e6) %>%
    ungroup() %>%
    inner_join(hla_annots, by = c("Name" = "tx_id")) %>%
    group_by(sampleid, gene_name) %>%
    summarise(cpm = sum(cpm)) %>%
    ungroup() %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_")

mapper <- "./hla-mapper/quant/%s_counts.txt" %>%
    sprintf(sample_ids) %>%
    setNames(sample_ids) %>%
    map_df(~read_tsv(., comment = "#"), .id = "sampleid") %>%
    select(sampleid = 1, tx_id = Geneid, counts = 7) %>%
    group_by(sampleid) %>%
    mutate(cpm = (counts/sum(counts)) * 1e6) %>%
    ungroup() %>%
    inner_join(hla_annots, by = "tx_id") %>%
    group_by(sampleid, gene_name) %>%
    summarise(cpm = sum(cpm)) %>%
    ungroup() %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_")

