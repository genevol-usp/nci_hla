devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(scales)
library(cowplot)
library(GGally)


# Data
samples <- readLines("./phase1/data/sample_ids.txt")

nci_qpcr <- read_tsv("./phase1/data/nci_expression.tsv")

nci_qpcr_gene <- distinct(nci_qpcr, subject, locus, mRNA, c_surface) %>%
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
    filter(subject %in% samples) 

nci_rnaseq <- 
    read_tsv("./phase1/3-map_to_transcriptome/hla_quantifications.tsv") %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    mutate(allele = hla_trimnames(gsub("IMGT_", "", allele), 3)) %>%
    select(subject, locus, allele, tpm)

nci_rnaseq_gene <- nci_rnaseq %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
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
  
nci_rnaseq_ref <- 
    read_tsv("./phase1/1-map_to_genome/hla_quantification.tsv") %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    select(subject, locus, tpm)

gene_df_ref <- inner_join(nci_rnaseq_ref, nci_qpcr_gene, by = c("subject", "locus"))


# plots
 
#Fig1 
png("./plots/Fig1.png", width = 8, height = 6, units = "cm", res = 300)
ggplot(nci_rnaseq_gene, aes(locus, tpm)) +
    geom_boxplot(fill = "grey") +
    theme_bw() +
    theme(text = element_text(size = 9, family = "Times"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = NULL, y = "TPM")
dev.off()


# Fig 2
png("./plots/Fig2.png", width = 10, height = 4, units = "cm", res = 300)
cor_df <- gene_df %>%
    group_by(locus) %>%
    do(data.frame(rho = cor(.$tpm, .$mRNA, method = "spearman"),
                  p = cor.test(.$tpm, .$mRNA, method = "spearman", exact = FALSE) %>% broom::tidy() %>% pull(p.value),
                  x = max(.$tpm),
                  y = max(.$mRNA))) %>%
    ungroup() %>%
    mutate(rho = round(rho, digits = 2),
           p = scientific(p, digits = 1),
           label = paste("rho ==", rho, "*', '~P ==", p)) 

ggplot(gene_df, aes(tpm, mRNA)) +
    geom_point(size = .75) +
    geom_smooth(method = lm, se = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    facet_wrap(~locus, scales = "free") +
    geom_label(data = cor_df, aes(x, y, label = label), 
               parse = TRUE, hjust = "inward", vjust = "inward", size = 2.5, 
               family = "Times",
               label.padding = unit(0.05, "lines"), 
               label.size = NA, alpha = 0.5) +
    theme_bw() +
    theme(text = element_text(size = 9, family = "Times"),
          strip.text = element_text(face = "bold")) +
    labs(y = "qPCR", x = "RNA-seq")
dev.off()

####

typing_errs <- read_tsv("./typings_nonconcordant.tsv") %>%
    distinct(subject, locus) 

gene_df_errs <- gene_df %>%
    anti_join(typing_errs)

cor_df_errs <- gene_df_errs %>%
    group_by(locus) %>%
    do(data.frame(rho = cor(.$tpm, .$mRNA, method = "spearman"),
                  p = cor.test(.$tpm, .$mRNA, method = "spearman", exact = FALSE) %>% broom::tidy() %>% pull(p.value),
                  x = max(.$tpm),
                  y = max(.$mRNA))) %>%
    ungroup() %>%
    mutate(rho = round(rho, digits = 2),
           p = scientific(p, digits = 1),
           label = paste("rho ==", rho, "*', '~P ==", p)) 

ggplot(gene_df_errs, aes(tpm, mRNA)) +
    geom_point(size = .75) +
    geom_smooth(method = lm, se = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    facet_wrap(~locus, scales = "free") +
    geom_label(data = cor_df_errs, aes(x, y, label = label), 
               parse = TRUE, hjust = "inward", vjust = "inward", size = 2.5, 
               family = "Times",
               label.padding = unit(0.05, "lines"), 
               label.size = NA, alpha = 0.5) +
    theme_bw() +
    theme(text = element_text(size = 9, family = "Times"),
          strip.text = element_text(face = "bold")) +
    labs(y = "qPCR", x = "RNA-seq")

####


# Fig S1
png("./plots/FigS1.png", width = 10, height = 4, units = "cm", res = 300)
cor_df_ref <- gene_df_ref %>%
    group_by(locus) %>%
    do(data.frame(rho = cor(.$tpm, .$mRNA, method = "spearman"),
                  p = cor.test(.$tpm, .$mRNA, method = "spearman", exact = FALSE) %>% broom::tidy() %>% pull(p.value),
                  x = max(.$tpm),
                  y = max(.$mRNA))) %>%
    ungroup() %>%
    mutate(rho = round(rho, digits = 2),
           p = scientific(p, digits = 1),
           label = paste("rho ==", rho, "*', '~P ==", p)) 

ggplot(gene_df_ref, aes(tpm, mRNA)) +
    geom_point(size = .75) +
    geom_smooth(method = lm, se = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    facet_wrap(~locus, scales = "free") +
    geom_label(data = cor_df_ref, aes(x, y, label = label), 
               parse = TRUE, hjust = "inward", vjust = "inward", size = 2.5, 
               family = "Times",
               label.padding = unit(0.05, "lines"), 
               label.size = NA, alpha = 0.5) +
    theme_bw() +
    theme(text = element_text(size = 9, family = "Times"),
          strip.text = element_text(face = "bold")) +
    labs(y = "qPCR", x = "RNA-seq")
dev.off()


# Fig 3
png("./plots/Fig3.png", width = 10, height = 4, units = "cm", res = 300)
cor_df_surface <- hlac_df %>%
    group_by(method) %>%
    do(data.frame(rho = cor(.$expression, .$c_surface, method = "spearman"),
                  p = cor.test(.$expression, .$c_surface, method = "spearman", exact = FALSE) %>% broom::tidy() %>% pull(p.value),
                  x = max(.$expression),
                  y = min(.$c_surface))) %>%
    ungroup() %>%
    mutate(rho = round(rho, digits = 2),
           p = scientific(p, digits = 1),
           label = paste("rho ==", rho, "*', '~P ==", p)) 


ggplot(hlac_df, aes(expression, c_surface)) +
    geom_smooth(method = lm, se = FALSE) +
    geom_point(size = .75) +
    facet_wrap(~method, scales = "free_x") +
    geom_label(data = cor_df_surface, aes(x, y, label = label), 
               parse = TRUE, hjust = "inward", vjust = "inward", size = 2.5, 
               family = "Times",
               label.padding = unit(0.05, "lines"), 
               label.size = NA, alpha = 0.5) +
    theme_bw() +
    labs(x = "RNA", y = "Surface Protein") +
    theme(text = element_text(size = 9, family = "Times"),
          strip.text = element_text(face = "bold")) 
dev.off()


# Fig 4
plot_lineages <- function(df) { 
    
    ggplot(data = df,
           aes(x = reorder(lineage, expression, FUN = median, na.rm = TRUE), 
               y = expression)) +
        ggbeeswarm::geom_quasirandom(aes(color = factor(hom)), 
                                     method = "smiley", varwidth = TRUE,
                                     size = .75, alpha = .75, show.legend = FALSE) +
        scale_color_manual(values = c("0" = "grey35", "1" = "red")) +
        geom_boxplot(outlier.shape = NA, fill = NA, alpha = .1) +
        scale_x_discrete(labels = function(x) sub("^(.+\\*)", "", x)) +
        scale_y_continuous(labels = comma, breaks = pretty_breaks(3)) +
        facet_wrap(~locus, scales = "free", ncol = 1, strip.position = "left") +
        labs(x = NULL, y = "mRNA") +
        theme_bw() +
        theme(text = element_text(size = 9, family = "Times"),
              strip.text = element_text(face = "bold"),
              legend.position = "none")
}

png("./plots/Fig4.png", width = 13.2, height = 8, units = "cm", res = 300)
p.rnaseq <- nci_rnaseq_lineage %>%
    filter(locus %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    group_by(lineage) %>%
    filter(n() >= 5) %>%
    ungroup() %>%
    plot_lineages() + labs(y = "TPM", title = "RNA-seq")

p.qpcr <- nci_qpcr_lineage %>%
    plot_lineages() + labs(y = "mRNA", title = "qPCR")

plot_grid(p.rnaseq, p.qpcr, ncol = 2)
dev.off()


#Fig 5
calc_trans_cors <- function(locus1, locus2, df) {
    
    m_sub <- df %>% 
        select(subject, locus1, locus2) %>%
        group_by(subject) %>%
        mutate_at(vars(locus2), rev) %>%
        ungroup() %>%
        select(-subject) %>%
        cor(use = "pairwise.complete.obs")
    
    m[c(locus1, locus2), c(locus1, locus2)] <<- m_sub 
}

hla_genes <- gencode_hla$gene_name

nci_rnaseq <- 
    read_tsv("./phase1/3-map_to_transcriptome/hla_quantifications.tsv") %>%
    mutate(allele = gsub("IMGT_", "", allele)) %>%
    select(subject, locus, allele, tpm) %>%
    filter(locus %in% hla_genes)

nci_rnaseq_gene <- nci_rnaseq %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup() 

global_cors <- nci_rnaseq_gene %>%
    mutate(locus = sub("HLA-", "", locus)) %>%
    spread(locus, tpm) %>% 
    select(!!! syms(sub("HLA-", "", gencode_hla$gene_name))) %>% 
    ggcorr(label = TRUE, label_round = 2, label_size = 2, hjust = 0.8, family = "Times") +
    labs(title = "Gene-level") +
    theme(legend.position = "none",
          text = element_text(size = 9, family = "Times"),
          plot.title = element_text(size = 9, family = "Times", hjust = 0.5))


phase <- "./phase1/4-phase_hla/phase_hla_haps.tsv" %>%
    read_tsv() %>%
    select(subject, locus, hap, allele, uncertain) %>%
    left_join(nci_rnaseq, by = c("subject", "locus", "allele")) %>%
    distinct() %>% 
    mutate(locus = sub("HLA-", "", locus)) %>%
    select(-uncertain, -allele) %>%
    spread(locus, tpm) %>%
    select(subject, !!! syms(sub("HLA-", "", gencode_hla$gene_name)))

cis_cors <- phase %>% 
    select(-subject) %>%
    ggcorr(label = TRUE, label_round = 2, label_size = 2, hjust = 0.8, family = "Times") +
    labs(title = "Within haplotypes") +
    theme(legend.position = "none",
          text = element_text(size = 9, family = "Times"),
          plot.title = element_text(size = 9, family = "Times", hjust = 0.5))

m <- 
    matrix(NA, nrow = length(hla_genes), ncol = length(hla_genes), 
           dimnames = list(sub("HLA-", "", hla_genes), 
                           sub("HLA-", "", hla_genes)))

phase_data <- tibble(locus1 = hla_genes, locus2 = hla_genes) %>% 
    mutate_all(~sub("HLA-", "", .)) %>% 
    expand(locus1, locus2) %>% 
    filter(locus1 != locus2) %>%
    pmap(calc_trans_cors, phase) 

trans_cors <- 
    ggcorr(data = NULL, cor_matrix = m, label = TRUE, label_round = 2,  
           label_size = 2, hjust = 0.8, family = "Times") +
    scale_fill_gradient2(name = "r", limits = c(0, 1), breaks = c(0, 0.5, 1), 
                         low = "#3B9AB2", mid = "#EEEEEE", high = "#F21A00") +
    labs(title = "Between haplotypes") +
    theme(text = element_text(size = 9, family = "Times"),
          plot.title = element_text(size = 9, family = "Times", hjust = 0.5))

legend <- get_legend(trans_cors)

trans_cors <- trans_cors + theme(legend.position = "none")

plot_correlations <- 
    plot_grid(global_cors, cis_cors, trans_cors, legend, nrow = 1, 
              rel_widths = c(1, 1, 1, .27))

classII_genes <-  c("HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1")

classII_and_CIITA <- gencode_chr_gene %>%
    filter(gene_name %in% c(classII_genes, "CIITA")) %>%
    arrange(as.integer(chr), start) %>%
    mutate(gene_name = factor(gene_name, levels = .$gene_name))

class_2_trans_df <- 
    "./phase1/3-map_to_transcriptome/gene_quantifications.tsv" %>%
    read_tsv() %>%
    inner_join(classII_and_CIITA, by = "gene_id") %>%
    select(subject, locus = gene_name, tpm) %>%
    spread(locus, tpm) %>%
    gather(locus, tpm, -subject, -CIITA) %>%
    mutate(locus = factor(locus, levels = classII_and_CIITA$gene_name)) %>%
    select(subject, locus, tpm, CIITA)

cor_df_classII <- class_2_trans_df %>%
    group_by(locus) %>%
    do(data.frame(rho = cor(.$tpm, .$CIITA, method = "spearman"),
                  x = min(.$tpm),
                  y = max(.$CIITA))) %>%
    ungroup() %>%
    mutate(rho = round(rho, digits = 2), 
           label = paste("rho ==", rho))

plot_ciita <- ggplot(class_2_trans_df, aes(tpm, CIITA)) +
    geom_point(size = .7) +
    geom_smooth(method = lm, se = FALSE) +
    scale_x_continuous(breaks = pretty_breaks(3)) +
    geom_label(data = cor_df_classII, aes(x, y, label = label), 
               parse = TRUE, hjust = "inward", vjust = "inward", size = 2.5, 
               family = "Times",
               label.padding = unit(0.05, "lines"), 
               label.size = NA, alpha = 0.5) +
    theme_bw() +
    theme(text = element_text(size = 9, family = "Times"),
          axis.text = element_text(hjust = 1)) +
    facet_wrap(~locus, nrow = 1, scales = "free") +
    labs(x = NULL)


png("./plots/Fig5.png", width = 13.2, height = 8, units = "cm", res = 300)
plot_grid(plot_correlations, plot_ciita, nrow = 2, ncol = 1, 
          rel_heights = c(1, .5), labels = c("A", "B"), label_size = 9, label_fontfamily = "Times",
          hjust = 0)
dev.off()


# Fig 6

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

png("./plots/Fig6.png", width = 12, height = 8, units = "cm", res = 300)
ggplot(summ_df, aes(q50.geuvadis, q50.nci, color = class)) +
    coord_cartesian(ylim = c(0, max(summ_df$q75.nci) + 50)) +
    ggrepel::geom_label_repel(aes(label = locus, fill = class), color = "white",
                              show.legend = FALSE, size = 2.5, 
                              label.padding = unit(0.05, "lines"), label.size = NA,
                              segment.color = "grey85") +
    geom_errorbar(aes(ymin = q25.nci, ymax = q75.nci)) +
    geom_errorbarh(aes(xmin = q25.geuvadis, xmax = q75.geuvadis, height = 30)) +
    geom_point(size = .75) +
    ggsci::scale_color_npg() +
    theme(text = element_text(size = 9, family = "Times")) +
    labs(x = "LCL", y = "PBL")
dev.off()



# Figs S2 and S3
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

png("./plots/FigS2.png", width = 10, height = 6, units = "cm", res = 300)
ggplot(timepoints_df, aes(tpm.t1, tpm.t2)) +
    geom_point(aes(color = locus)) +
    ggsci::scale_color_npg() +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    theme_bw() +
    theme(text = element_text(size = 9, family = "Times")) +
    labs(x = "TPM (timepoint 1)", y = "TPM (timepoint 2)")
dev.off()


png("./plots/FigS3.png", width = 10, height = 10, units = "cm", res = 300)
cor_df_timepoints <- timepoints_df %>%
    group_by(locus) %>%
    do(data.frame(rho = cor(.$tpm.t1, .$tpm.t2, method = "spearman"),
                  x = max(.$tpm.t1),
                  y = min(.$tpm.t2))) %>%
    ungroup() %>%
    mutate(rho = round(rho, digits = 2), 
           label = paste("rho ==", rho))


ggplot(timepoints_df, aes(tpm.t1, tpm.t2)) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_point() +
    geom_label(data = cor_df_timepoints, aes(x, y, label = label), 
               parse = TRUE, hjust = "inward", vjust = "inward", size = 2.5, 
               family = "Times",
               label.padding = unit(0.05, "lines"), 
               label.size = NA, alpha = 0.5) +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    facet_wrap(~locus, scales = "free") +
    theme_bw() +
    theme(text = element_text(size = 9, family = "Times")) +
    labs(x = "timepoint 1", y = "timepoint 2")
dev.off()









# B2M-normalized
b2m <- read_tsv("./phase1/3-map_to_transcriptome/b2m_expression.tsv") %>%
    select(subject, b2m_counts = counts)

qpcr <- nci_qpcr_gene %>% select(subject, locus, qPCR = mRNA)

b2m_norm_df <-
    read_tsv("./phase1/3-map_to_transcriptome/hla_quantifications.tsv") %>%
    filter(locus %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    group_by(subject, locus) %>%
    summarise(est_counts = sum(est_counts)) %>%
    ungroup() %>%
    left_join(b2m, by = "subject") %>%
    mutate(norm_counts = est_counts/b2m_counts) %>%
    select(subject, locus, norm_counts) %>%
    left_join(qpcr, by = c("subject", "locus")) %>%
    filter(norm_counts < 10)

png("./plots/b2mnorm_vs_qpcr.png", width = 8, height = 3, units = "in", res = 300)
ggplot(b2m_norm_df, aes(qPCR, norm_counts)) +
    geom_point(size = 2) +
    geom_smooth(method = lm, se = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    facet_wrap(~locus, scales = "free") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          strip.text = element_text(face = "bold", size = 16)) +
    labs(x = "qPCR", y = "RNAseq (B2M-normalized counts)") +
    ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 2,
                          formula = y ~ x, parse = TRUE, size = 6)
dev.off()


####









# Timepoint 2 vs qPCR
gene_df_t2 <- 
    inner_join(nci_rnaseq_gene_t2, nci_qpcr_gene, by = c("subject", "locus"))

png("./plots/rnaseqt2_vs_qpcr.png", width = 12, height = 4, units = "in", res = 300)
ggplot(gene_df_t2, aes(mRNA, tpm)) +
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









