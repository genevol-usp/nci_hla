library(tidyverse)
library(tidytext)
library(scales)
library(cowplot)
library(ggsci)

# read length and coverage
ids_reads_t1 <- dir("/raid/genevol/nci_rnaseq/phase1/read_info") %>%
    sub("\\.txt$", "", .)

reads_t1 <- sprintf("/raid/genevol/nci_rnaseq/phase1/read_info/%s.txt", ids_reads_t1) %>%
    setNames(ids_reads_t1) %>%
    map_df(read_lines, .id = "sampleid") %>%
    pivot_longer(everything(), names_to = "sampleid") %>%
    separate(value, c("depth", "readlen"), sep = " ", convert = TRUE) %>%
    mutate(sampleid = str_extract(sampleid, "(66K\\d+)")) %>%
    group_by(sampleid) %>%
    summarise(depth = sum(depth),
              readlen = mean(readlen)) %>%
    ungroup()

ids_reads_t2 <- dir("/raid/genevol/nci_rnaseq/phase2/read_info") %>%
    sub("\\.txt$", "", .)

reads_t2 <- sprintf("/raid/genevol/nci_rnaseq/phase2/read_info/%s.txt", ids_reads_t2) %>%
    setNames(ids_reads_t2) %>%
    map_df(read_lines, .id = "sampleid") %>%
    pivot_longer(everything(), names_to = "sampleid") %>%
    separate(value, c("depth", "readlen"), sep = " ", convert = TRUE) %>%
    group_by(sampleid) %>%
    summarise(depth = sum(depth),
              readlen = mean(readlen)) %>%
    ungroup()

reads_df <- bind_rows("Timepoint 1" = reads_t1,
                      "Timepoint 2" = reads_t2, .id = "timepoint")

p1 <- ggplot(reads_df, aes(reorder_within(sampleid, depth, timepoint), depth)) +
    geom_col(width = 1) +
    scale_y_continuous(labels = function(x) round(floor(x)/1e6),
                       breaks = pretty_breaks(10)) +
    facet_grid(~timepoint, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()) +
    labs(x = NULL, y = NULL)

p2 <- ggplot(reads_df, aes(sampleid, readlen)) +
    geom_col(width = 1) +
    facet_grid(~timepoint, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()) +
    labs(x = NULL, y = NULL,
         title = "Read length")



tx_annots <- read_tsv("../indices/transcript_annotation_df.tsv") %>%
    select(tx_id, gene_id, gene_name, tx_type)

ids <- dir("./salmon/quant")

salmon_tx <- "./salmon/quant/%s/quant.sf" %>%
    sprintf(ids) %>%
    setNames(ids) %>%
    map_df(. %>% read_tsv %>% select(tx_id = Name, counts = NumReads, tpm = TPM), 
           .id = "sampleid") %>%
    mutate(timepoint = ifelse(grepl("t1$", sampleid), "Timepoint 1", "Timepoint 2")) %>%
    left_join(tx_annots) %>%
    select(sampleid, timepoint, tx_id, tx_type, gene_id, gene_name, counts, tpm)

meantpm_by_types <- salmon_tx %>%
    group_by(sampleid, timepoint, tx_type) %>%
    summarise(tpm = sum(tpm)) %>%
    group_by(sampleid, timepoint) %>%
    mutate(tpm = tpm/sum(tpm)) %>%
    group_by(timepoint, tx_type) %>%
    summarise(prop = mean(tpm)) %>%
    group_by(tx_type) %>%
    mutate(prop_type = sum(prop)) %>%
    ungroup() %>%
    mutate(tx_type = ifelse(prop_type > 0.01, tx_type, "Other")) %>%
    group_by(timepoint, tx_type) %>%
    summarise(prop = mean(prop)) %>%
    ungroup() %>%
    arrange(desc(prop)) %>%
    mutate(tx_type = recode(tx_type, "nonsense_mediated_decay" = "NMD"))

type_order <- meantpm_by_types %>%
    filter(timepoint == "Timepoint 1") %>%
    arrange(prop) %>%
    pull(tx_type)

cols <- c("grey", pal_npg()(10), "black")

p3 <- meantpm_by_types %>%
    mutate(tx_type = factor(tx_type, levels = type_order)) %>%
    ggplot(aes(timepoint, prop, fill = tx_type)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = cols) +
    scale_y_continuous(labels = percent) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    guides(fill = guide_legend(keyheight = .5)) +
    labs(x = NULL, y = NULL,
         fill = "Biotype")



salmon_gene <- salmon_tx %>%
    group_by(sampleid, timepoint, gene_id, gene_name) %>%
    summarise(tpm = sum(tpm)) %>%
    ungroup()

tpm_proportions <- salmon_gene %>%
    group_by(sampleid, timepoint) %>%
    mutate(prop = tpm/sum(tpm)) %>%
    ungroup() %>%
    mutate(labl = ifelse(prop >= 0.05, gene_name, "Others")) %>%
    group_by(sampleid, timepoint, gene = labl) %>%
    summarise(prop = sum(prop)) %>%
    ungroup()

gene_levels <- tpm_proportions %>%
    group_by(gene) %>%
    summarise(prop = mean(prop)) %>%
    ungroup() %>%
    arrange(-prop) %>%
    pull(gene)

others_prop <- tpm_proportions %>%
    filter(gene == "Others") %>%
    arrange(prop)

tpm_proportions <- tpm_proportions %>%
    mutate(sampleid = factor(sampleid, levels = others_prop$sampleid),
           gene = factor(gene, levels = rev(gene_levels)),
           timepoint = recode(timepoint, "t1" = "Timepoint 1", "t2" = "Timepoint 2"))


mycols <- rev(c("grey", pal_npg()(9)))

p4 <- ggplot(tpm_proportions, aes(sampleid, prop, fill = gene)) +
    geom_bar(stat = "identity", position = "fill", width = 1) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = mycols) +
    facet_grid(~timepoint, scales = "free_x", space = "free_x", switch = "x") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank()) +
    guides(fill = guide_legend(override.aes = list(color = "black"),
                               reverse = TRUE, 
                               keyheight = .7)) +
    labs(y = NULL, fill = "Gene")

# comparison with Geuvadis
salmon_tx %>%
    group_by(sampleid, timepoint) %>%
    mutate(prop = counts/sum(counts)) %>%
    ungroup() %>% 
    filter(gene_name %in% c("RN7SL1", "RN7SL2", "RN7SK")) %>%
    group_by(timepoint, gene_id, gene_name) %>%
    summarise(prop = mean(prop)) %>%
    group_by(timepoint) %>%
    summarise(sum(prop)) %>%
    ungroup()
    

salmon_gene %>%
    group_by(sampleid, timepoint) %>%
    mutate(prop = tpm/sum(tpm)) %>%
    ungroup()  %>% 
    filter(gene_name %in% c("RN7SL1", "RN7SL2", "RN7SK")) %>%
    group_by(timepoint, gene_id, gene_name) %>%
    summarise(prop = mean(prop)) %>%
    group_by(timepoint) %>%
    summarise(sum(prop)) %>%
    ungroup()

library(furrr)

plan(multisession, workers = 16)

geuvadis_ids <- dir("/raid/genevol/geuvadis/salmon/quant/")

salmon_geuv <- "/raid/genevol/geuvadis/salmon/quant/%s/quant.sf" %>%
    sprintf(geuvadis_ids) %>%
    setNames(geuvadis_ids) %>%
    future_map_dfr(read_tsv, .id = "sampleid") %>%
    select(sampleid, tx_id = Name, tpm = TPM, counts = NumReads)

salmon_geuv %>%
    group_by(sampleid) %>%
    mutate(prop = counts/sum(counts)) %>%
    ungroup() %>% 
    filter(tx_id %in% txs) %>%
    group_by(tx_id) %>%
    summarise(prop = mean(prop)) %>%
    ungroup()

###########################


# PCA on phenotypes
sample_info <- "/raid/genevol/nci_rnaseq/phase1/nci_sample_info.tsv" %>%
    read_tsv() %>%
    select(sampleid = sample_id, batch) %>%
    distinct()

depth_df <- reads_df %>%
    filter(timepoint == "Timepoint 1") %>%
    select(sampleid, depth)

gene_props <- salmon_gene %>%
    filter(timepoint == "Timepoint 1") %>%
    mutate(sampleid = str_remove(sampleid, "_t1")) %>%
    group_by(sampleid, timepoint) %>%
    mutate(prop = tpm/sum(tpm)) %>%
    ungroup()

top_genes <- gene_props %>%
    group_by(gene_id, gene_name) %>%
    summarise(m = mean(prop)) %>%
    ungroup() %>%
    arrange(desc(m)) %>%
    slice(1:3) %>%
    pull(gene_id)

rnp_df <- gene_props %>%
    filter(gene_id %in% top_genes) %>%
    group_by(sampleid) %>%
    summarise(rnp = sum(prop)) %>%
    ungroup()


pca_pheno <- read_delim("./salmon/out.pca", delim = " ") %>%
    pivot_longer(-1, names_to = "sampleid") %>%
    rename(pc = SampleID) %>%
    mutate(pc = str_remove(pc, "out_1_1_svd_")) %>%
    pivot_wider(names_from = pc, values_from = value) %>%
    left_join(sample_info, by = "sampleid") %>%
    left_join(depth_df, by = "sampleid") %>%
    left_join(others_df, by = "sampleid") %>%
    select(sampleid, batch, depth, rnp, everything())

p_pca_1 <- ggplot(pca_pheno, aes(PC1, PC2)) +
    geom_point(aes(color = rnp), size = 2) +
    scale_color_viridis_c() +
    theme(panel.background = element_rect(fill = "grey50"),
          panel.grid = element_blank(),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10, margin = unit(c(0, -.1, 0, 0), "cm")),
          axis.text.x = element_text(hjust = 1)) +
    labs(color = "RNP RNAs")

p_pca_2 <- ggplot(pca_pheno, aes(PC2, PC3)) +
    geom_point(aes(color = rnp), size = 2) +
    scale_color_viridis_c() +
    guides(color = guide_colorbar(barwidth = .5)) +
    theme(panel.background = element_rect(fill = "grey50"),
          panel.grid = element_blank(),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10, margin = unit(c(0, -.1, 0, 0), "cm")),
          axis.text.x = element_text(hjust = 1)) +
    labs(color = "RNP\nRNAs")

p5 <- plot_grid(
    p_pca_1 + theme(legend.position = "none"), 
    p_pca_2 + theme(legend.position = "none"),
    get_legend(p_pca_2), 
    NULL,
    nrow = 1, rel_widths = c(1, 1, .2, .1), 
    greedy = TRUE
)

plot_grid(
    NULL,
    plot_grid(p1, NULL, nrow = 1, rel_widths = c(1, .2)),
    NULL,
    p4,
    NULL,
    plot_grid(NULL, p3, NULL, p5, nrow = 1,
              labels = c("C)", "", "D)", NULL), 
              rel_widths = c(.075, .9, .075, 1)),
    ncol = 1, 
    rel_heights = c(.1, .75, .1, 1, .1, .6),
    labels = c("A)", "", "B)", "", NULL, NULL)
)

ggsave("./plots/diagnostics.jpeg", width = 8.5, height = 8)
