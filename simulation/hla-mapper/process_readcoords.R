library(tidyverse)

# annotations
hla_txs <- "../make_simulation_data/simulation_index.fa" %>%
    Biostrings::readDNAStringSet() %>%
    names() %>%
    grep("_[ABC]\\*", ., value = TRUE) %>%
    str_extract("(ENST[0-9.]+)") %>%
    unique()

annots <- "~/gencode/gencode.v37.primary_assembly.annotation.gtf" %>%
    read_tsv(comment = "#", col_names = FALSE, col_types = "c-cii-c-c")

exon_annots <- annots %>%
    filter(X1 == "chr6", X3 == "exon", X4 > 29e6, X5 < 33e6) %>%
    transmute(start = X4, end = X5, info = X9, 
	      tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+")) %>%
    filter(tx_id %in% hla_txs) %>%
    transmute(start, end, tx_id, 
	      gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+"),
	      exon_number = str_extract(info, "(?<=exon_number\\s)\\d+")) %>%
    mutate(pos = map2(start, end, ~.x:.y)) %>%
    select(-start, -end) %>%
    unnest(c(pos)) %>%
    group_by(tx_id) %>%
    mutate(i = 1:n()) %>%
    ungroup()


# sample IDs
samples <- readLines("../samples.txt")

# reads
readmaps <- read_rds("./simulated_reads_mappings.rds") %>%
    separate(readname, c("readid", "tx_id", "allele"), sep = "_") %>%
    select(sampleid, readid, tx_id, mapped_gene)

readcoords <- sprintf("~/simulation/readcoords/%s.txt", samples) %>%
    setNames(samples)  %>%
    map_df(~read_tsv(., col_names = FALSE), .id = "sampleid") %>%
    extract(X1, c("readid", "tx_id", "allele", "coords"), "(read\\d+)_(ENST[0-9.]+)_(\\S+) (.+)") %>%
    separate(coords, c("mate1", "mate2"), sep = ";") %>%
    mutate(across(mate1:mate2, ~sub("^mate[12]:(\\d+-\\d+).*$", "\\1", .)))


summarize_mappings <- function(sample_id) {

    readcoords_i <- readcoords %>% 
	filter(sampleid == sample_id) %>%
	select(-allele)

    readmaps_i <- readmaps %>% 
	filter(sampleid == sample_id)

    tmp <- left_join(readcoords_i, readmaps_i, by = c("sampleid", "readid", "tx_id")) %>%
	pivot_longer(mate1:mate2, names_to = "mate") %>%
	separate(value, c("start", "end"), sep = "-", convert = TRUE) %>%
	mutate(i = map2(start, end, ~.x:.y)) %>%
	select(-start, -end) %>%
	unnest(c(i)) %>%
	left_join(distinct(exon_annots, tx_id, gene_name), by = "tx_id")

    out <- tmp %>% 
	count(sampleid, gene_name, mapped_gene, i) %>%
	arrange(sampleid, gene_name, i, mapped_gene) %>%
	group_by(sampleid, gene_name, i) %>%
	mutate(prop = n/sum(n)) %>%
	ungroup()
	
    out
}

library(furrr)

plan(multisession, workers = 10)

results <- future_map_dfr(samples, summarize_mappings)

results_baselevel <- results %>%
    group_by(gene_name, mapped_gene, i) %>%
    summarise(prop = mean(prop)) %>%
    ungroup() %>%
    mutate(mapped_gene = replace_na(mapped_gene, "Not present in SAM")) %>%
    arrange(gene_name, i, desc(prop))

write_rds(results_baselevel, "../plot_data/summary_subread_mappings.rds")

