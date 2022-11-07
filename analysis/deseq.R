library(tidyverse)
library(tximport)
library(DESeq2)

fix_quants <- function(f) {

    fin <- read_tsv(f)
    
    non_hla <- filter(fin, !grepl("_[ABC]\\*", Name))
    hla <- filter(fin, grepl("_[ABC]\\*", Name)) 

    hla_summ <- hla %>%
	mutate(Name = sub("^([^_]+).*$", "\\1", Name)) %>%
	group_by(Name) %>%
	summarise(Length = ifelse(any(NumReads > 0), weighted.mean(Length, NumReads), mean(Length)),
		  EffectiveLength = ifelse(any(NumReads > 0), weighted.mean(EffectiveLength, NumReads), mean(EffectiveLength)),
		  TPM = sum(TPM),
		  NumReads = sum(NumReads)) %>%
	ungroup()

    bind_rows(non_hla, hla_summ)
}

samples_df <- read_tsv("./samples.txt", col_names = c("batch", "sampleid"))

samples_tps <- samples_df %>%
    group_by(sampleid) %>%
    filter(all(1:2 %in% batch)) %>%
    ungroup() %>%
    unite("id", c("sampleid", "batch"), sep = "_t") %>%
    pull(id)

files <- file.path("./salmon-pers/quant", samples_tps, "quant.sf") %>%
    setNames(samples_tps)

fix_quants(files[1]) %>% tail

out_dir <- "./salmon-pers/quant_hlasumm" 
dir.create(out_dir)
walk(samples_tps, ~file.path(out_dir, .x) %>% dir.create())

walk(samples_tps, ~fix_quants(files[.x]) %>% 
     write_tsv(file.path("./salmon-pers/quant_hlasumm", .x, "quant.sf")))

files_salmon <- file.path("./salmon-pers/quant_hlasumm", samples_tps, "quant.sf") %>%
    setNames(samples_tps)

tx_annots <- read_tsv("../indices/transcript_annotation_df.tsv")

gene_annots <- distinct(tx_annots, gene_id, gene_name)

tx_to_gene <- tx_annots %>%
    select(TXNAME = tx_id, GENEID = gene_id)

transcs <- read_tsv(files_salmon[1]) %>% select(TXNAME = Name)
transcs2 <- read_tsv(files_salmon[11]) %>% select(TXNAME = Name)

tx_to_gene_adj <- inner_join(transcs, tx_to_gene)

txi <- tximport(files_salmon, type = "salmon", tx2gene = tx_to_gene_adj)

col_data <- tibble(sampleid = samples_tps) %>%
    separate(sampleid, c("donor", "batch"), sep = "_") %>%
    mutate_all(factor)

dds <- DESeqDataSetFromTximport(txi, col_data, ~donor + batch)
dds <- estimateSizeFactors(dds)
ncts <- counts(dds, normalized=TRUE)

out <- as_tibble(ncts, rownames = "gene_id") %>%
    left_join(gene_annots, by = "gene_id") %>%
    select(gene_id, gene_name, everything()) %>%
    pivot_longer(starts_with("66K"), names_to = "sampleid", values_to = "scaled_counts") %>%
    select(sampleid, gene_id, gene_name, scaled_counts) %>%
    group_by(sampleid) %>%
    nest() %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_") %>%
    unnest(cols = c(data))

write_tsv(out, "./plot_data/deseq_counts.tsv")


# normalize separate
txi_1 <- tximport(files_salmon[1:11], type = "salmon", tx2gene = tx_to_gene_adj)

col_data_1 <- tibble(x = 1:11)

dds_1 <- DESeqDataSetFromTximport(txi_1, col_data_1, ~1)
dds_1 <- estimateSizeFactors(dds_1)
ncts_1 <- counts(dds_1, normalized = TRUE)

out_1 <- as_tibble(ncts_1, rownames = "gene_id") %>%
    left_join(gene_annots, by = "gene_id") %>%
    select(gene_id, gene_name, everything()) %>%
    pivot_longer(starts_with("66K"), names_to = "sampleid", values_to = "scaled_counts") %>%
    select(sampleid, gene_id, gene_name, scaled_counts) %>%
    group_by(sampleid) %>%
    nest() %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_") %>%
    unnest(cols = c(data))

txi_2 <- tximport(files_salmon[12:22], type = "salmon", tx2gene = tx_to_gene_adj)

col_data_2 <- tibble(x = 1:11)

dds_2 <- DESeqDataSetFromTximport(txi_2, col_data_2, ~1)
dds_2 <- estimateSizeFactors(dds_2)
ncts_2 <- counts(dds_2, normalized = TRUE)

out_2 <- as_tibble(ncts_2, rownames = "gene_id") %>%
    left_join(gene_annots, by = "gene_id") %>%
    select(gene_id, gene_name, everything()) %>%
    pivot_longer(starts_with("66K"), names_to = "sampleid", values_to = "scaled_counts") %>%
    select(sampleid, gene_id, gene_name, scaled_counts) %>%
    group_by(sampleid) %>%
    nest() %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_") %>%
    unnest(cols = c(data))

bind_rows(out_1, out_2) %>%
    write_tsv("./plot_data/deseq_counts_separatenorm.tsv")





write_tsv(out, "./plot_data/deseq_counts.tsv")





