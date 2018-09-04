devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(scales)
library(cowplot)

# Functions
plot_lineages <- function(data) {
    ggplot(data, aes(x = reorder(lineage, expression, FUN = mean, na.rm = TRUE), 
		     y = expression)) +
	ggbeeswarm::geom_quasirandom(aes(color = factor(hom)), 
				     method = "smiley", varwidth = TRUE,
				     size = .75, alpha = .75, show.legend = FALSE) +
    scale_x_discrete(labels = function(x) sub("^(.+\\*)", "", x)) +
    scale_y_continuous(labels = comma, breaks = pretty_breaks(3)) +
	scale_color_manual(values = c("0" = "grey25", "1" = "green")) +
	stat_summary(fun.y = mean, geom = "point", shape = "\U2014", size = 9) +
	facet_wrap(~locus, scales = "free", ncol = 1, strip.position = "left") +
	labs(x = NULL, y = "mRNA") +
	theme_bw() +
	theme(axis.title = element_text(size = 16),
	      axis.text.y = element_text(size = 14),
	      axis.text.x = element_text(size = 11),
	      legend.text = element_text(size = 12),
	      strip.text = element_text(face = "bold", size = 14))
}

# Data
samples <- readLines("./phase1/data/sample_ids.txt")

nci_qpcr <- read_tsv("./phase1/data/nci_expression.tsv")

nci_qpcr_gene <- distinct(nci_qpcr_ab, subject, locus, mRNA, c_surface) %>%
    filter(subject %in% samples)

nci_qpcr_lineage <- nci_qpcr %>%
    mutate(lineage = hla_trimnames(allele, 1)) %>%
    select(subject, locus, lineage, expression = mRNA) %>%
    group_by(subject, locus) %>%
    mutate(hom = ifelse(length(unique(lineage)) == 1, 1L, 0L)) %>%
    filter(!is.na(expression)) %>%
    group_by(lineage) %>%
    filter(n() > 5) %>%
    ungroup() %>%
    mutate(est_allele_exp = lm(expression ~ lineage - 1)$fitted.values) %>%
    filter(subject %in% samples) %>%
    group_by(lineage) %>%
    filter(n() > 5) %>%
    ungroup()

nci_rnaseq <- 
    read_tsv("./phase1/3-map_to_transcriptome/hla_quantifications.tsv") %>%
    mutate(allele = hla_trimnames(gsub("IMGT_", "", allele), 3)) %>%
    select(subject, locus, allele, tpm)

nci_rnaseq_gene <- nci_rnaseq %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup() %>%
    group_by(locus) %>%
    filter(mean(tpm) >= 10) %>%
    ungroup() %>%
    mutate(locus = reorder(locus, tpm, median),
           locus = factor(locus, levels = rev(levels(locus))))

nci_rnaseq_lineage <- nci_rnaseq %>%
    mutate(lineage = hla_trimnames(allele, 1)) %>%
    group_by(subject, locus) %>%
    mutate(hom = ifelse(length(unique(allele)) == 1, 1L, 0L)) %>%
    select(subject, locus, lineage, hom, expression = tpm) %>%
    ungroup()

gene_df <- inner_join(nci_rnaseq_gene, nci_qpcr_gene, by = c("subject", "locus"))

hlac_df <- gene_df %>%
    filter(locus == "HLA-C", !is.na(c_surface)) %>%
    select(subject, rnaseq = tpm, qPCR = mRNA, c_surface) %>%
    gather(method, expression, rnaseq, qPCR) %>%
    select(subject, method, expression, c_surface) %>%
    arrange(subject, method)
  

# plots
png("./plots/expression_boxplot.png", width = 8, height = 5, units = "in", res = 300)
ggplot(nci_rnaseq_gene, aes(locus, tpm)) +
    geom_boxplot(fill = "grey") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 12)) +
    labs(y = "TPM")
dev.off()

png("./plots/rnaseq_vs_qpcr.png", width = 12, height = 4, units = "in", res = 300)
ggplot(gene_df, aes(mRNA, tpm)) +
    geom_point(size = 2) +
    geom_smooth(method = lm, se = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    facet_wrap(~locus, scales = "free") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          strip.text = element_text(face = "bold", size = 16)) +
    labs(x = "qPCR", y = "RNAseq (TPM)") +
    ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 2,
                          formula = y ~ x, parse = TRUE, size = 6)
dev.off()

png("./plots/ab_vs_rna.png", width = 10, height = 4, units = "in", res = 300)
ggplot(hlac_df, aes(expression, c_surface)) +
  geom_point(size = 2) +
  facet_wrap(~method, scales = "free") +
  theme_bw() +
  labs(x = "RNA", y = "Surface Protein") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(face = "bold", size = 16)) +
    ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 2,
			  formula = y ~ x, parse = TRUE, size = 6)
dev.off()

png("./plots/lineages.png", width = 10, height = 8, units = "in", res = 300)
p.rnaseq <- nci_rnaseq_lineage %>%
    filter(locus %in% paste0("HLA-", c("A", "B", "C", "DPB1", "DQB1", "DRB1"))) %>%
    group_by(lineage) %>%
    filter(n() > 5) %>%
    ungroup() %>%
    plot_lineages() + labs(y = "TPM", title = "RNA-seq")

p.qpcr <- nci_qpcr_lineage %>%
    plot_lineages() + labs(y = "mRNA (qPCR)", title = "qPCR")

left_plot <- plot_grid(p.qpcr, NULL, nrow = 2, rel_heights = c(1, .9))

plot_grid(p.rnaseq, left_plot, ncol = 2)
dev.off()

nci_rnaseq_t2 <- 
    read_tsv("./phase2/3-map_to_transcriptome/hla_quantifications.tsv") %>%
    mutate(allele = hla_trimnames(gsub("IMGT_", "", allele), 3)) %>%
    select(subject, locus, allele, tpm)

nci_rnaseq_gene_t2 <- nci_rnaseq_t2 %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup() %>%
    mutate(locus = reorder(locus, tpm, median),
           locus = factor(locus, levels = rev(levels(locus))))

timepoints_df <- 
    inner_join(nci_rnaseq_gene, nci_rnaseq_gene_t2, 
               by = c("subject", "locus"), suffix = c(".t1", ".t2"))

png("./plots/nci_timepoints.png", width = 6, height = 4, units = "in", res = 300)
ggplot(timepoints_df, aes(tpm.t1, tpm.t2)) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_point() +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    ggpmisc::stat_poly_eq(formula = y ~ x, 
                          aes(label = ..adj.rr.label..), parse = TRUE, size = 3.5,
                          label.x.npc = .8, label.y.npc = 0, rr.digits = 2) +
    facet_wrap(~locus, scales = "free") +
    theme_bw() +
    labs(x = "timepoint 1", y = "timepoint 2")
dev.off()

# Comparison with Geuvadis

gene_nci <- nci_rnaseq_gene %>%
    mutate(locus = as.character(locus)) %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    split(.$locus) %>%
    map_df(. %>% pull(tpm) %>% fivenum() %>% .[2:4] %>% t() %>% as.data.frame(),
           .id = "locus") %>%
    rename(q25 = V1, q50 = V2, q75 = V3)

gene_geuvadis <- 
    "../hlaexpression/geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    select(subject, locus, allele, tpm) %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    group_by(subject, locus) %>%
    summarise(tpm = sum(tpm)) %>%
    ungroup() %>%
    split(.$locus) %>%
    map_df(. %>% pull(tpm) %>% fivenum() %>% .[2:4] %>% t() %>% as.data.frame(),
           .id = "locus") %>%
    rename(q25 = V1, q50 = V2, q75 = V3)

summ_df <- 
    inner_join(gene_geuvadis, gene_nci, by = "locus", 
               suffix = c(".geuvadis", ".nci")) %>%
    mutate(locus = sub("HLA-", "", locus),
           class = ifelse(locus %in% c("A", "B", "C"), "Class I", "Class II"),
           ny = ifelse(class == "Class I", max(q75.nci) + 50 - q50.nci, 0 - q50.nci))

png("./plots/nci_vs_geuvadis_summaries.png", width = 6, height = 4, units = "in", res = 150)
ggplot(summ_df, aes(q50.geuvadis, q50.nci, color = class)) +
    coord_cartesian(ylim = c(0, max(summ_df$q75.nci) + 50)) +
    ggrepel::geom_label_repel(aes(label = locus, fill = class), color = "white",
                              show.legend = FALSE, size = 3.5, 
                              nudge_y = summ_df$ny,
                              segment.color = "grey75") +
    geom_errorbar(aes(ymin = q25.nci, ymax = q75.nci)) +
    geom_errorbarh(aes(xmin = q25.geuvadis, xmax = q75.geuvadis, height = 30)) +
    geom_point(size = 3) +
    ggsci::scale_color_npg() +
    labs(x = "LCL", y = "PBMC", title = "Median, 25th and 75th quantiles of expression")
dev.off()

