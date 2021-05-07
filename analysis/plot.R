library(tidyverse)
library(scales)
library(ggthemes)
library(ggbeeswarm)
library(cowplot)
library(tidytext)
library(ggrepel)


salmon <- read_rds("./plot_data/salmon_quants.rds")

tpm_proportions <- salmon %>%
    arrange(sampleid, -tpm) %>%
    group_by(sampleid) %>%
    mutate(prop = tpm/sum(tpm)) %>%
    ungroup() %>%
    mutate(labl = ifelse(prop >= 0.05, gene_name, "Others")) %>%
    group_by(sampleid, gene = labl) %>%
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
           gene = factor(gene, levels = gene_levels))

mycols <- c("grey", "#046C9A", "#D69C4E", "#ABDDDE",
            "tomato3", "darkgreen", "grey35", "#ECCBAE")

p1 <- ggplot(tpm_proportions, aes(sampleid, prop, fill = gene)) +
    geom_bar(stat = "identity", position = "fill", width = 0.95, alpha = .75,
             color = "grey85") +
    geom_text(aes(label = sub("%", "", scales::percent(prop, accuracy = 1))),
              size = 2.5, color = "grey25", fontface = "bold",
              position = position_stack(vjust = .85)) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = mycols,
		      guide = guide_legend(direction = "horizontal")) +
    coord_flip() +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank()) +
    guides(fill = guide_legend(override.aes = list(color = "black"))) +
    labs(y = "% of whole gene expression", fill = "Gene")
          

ggsave("./plots/expression_proportions.jpeg", p1, width = 8, height = 10)


### Salmon pers vs qPCR
#### locus overall expression
salmon_pers_tpm <- read_tsv("./salmon-pers/quants.bed.gz") 

salmon_pers_hla <- salmon_pers_tpm %>%
    filter(id %in% c("HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", 
                     "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1")) %>%
    pivot_longer(-gene_name, names_to = "sampleid", values_to = "rna") %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_") %>%
    filter(timepoint == "t1") %>%
    select(sampleid, gene_name, rna)
    
ggplot(salmon_pers_tpm, aes(fct_reorder(gene_name, rna, .desc = TRUE), rna)) +
    geom_boxplot(fill = "grey") +
    theme_bw() +
    labs(x = NULL, y = "TPM")

ggsave("./plots/salmon_pers_gene.jpg", width = 6.5, height = 2.5)


#### Salmon vs qPCR

nci_expression <- 
    "/raid/genevol/nci_rnaseq/phase1/HLA-A, -B, -C expression levels for Pat.xlsx" %>%
    readxl::read_excel() %>%
    janitor::clean_names() %>%
    select(sampleid = 3, ends_with("rna"), ends_with("surface")) %>%
    mutate(across(-sampleid, as.numeric)) %>%
    pivot_longer(-sampleid, names_to = "tmp") %>%
    extract(tmp, c("gene_name", "unit"), "(hla_[abc])_(.+)") %>%
    mutate(gene_name = toupper(gene_name),
           gene_name = sub("_", "-", gene_name)) %>%
    pivot_wider(names_from = unit, values_from = value)
    
    
salmon_pers <- read_tsv("./salmon-pers/quants_std.bed") %>%
    filter(id %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    select(gene_name = id, starts_with("66K")) %>%
    pivot_longer(-gene_name, names_to = "sampleid", values_to = "rnaseq") %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_")

quant_df <- salmon_pers %>%
    filter(timepoint == "t1") %>%
    select(-timepoint) %>%
    left_join(nci_expression, by = c("sampleid", "gene_name")) %>%
    select(sampleid, gene_name, qPCR = m_rna, rnaseq, surface)

cor_df <- quant_df %>%
    group_by(gene_name) %>%
    summarise(y = max(qPCR),
              r = round(cor(qPCR, rnaseq), 2),
              rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

ggplot(quant_df, aes(rnaseq, qPCR)) +
    geom_point(size = .75) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_text(data = cor_df, aes(x = -2.5, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3) +
    facet_wrap(~gene_name, scales = "free_y") +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "RNA-seq")

ggsave("./plots/real_data.jpeg", height = 2, width = 6)

# Salmon ref and hla-mapper vs qPCR

## Salmon
salmon_ref <- read_tsv("./salmon/quants_std.bed") %>%
    filter(id %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    select(gene_name = id, starts_with("66K")) %>%
    pivot_longer(-gene_name, names_to = "sampleid", values_to = "rnaseq")

quant_ref_df <- salmon_ref %>%
    left_join(nci_expression, by = c("sampleid", "gene_name")) %>%
    select(sampleid, gene_name, qPCR = m_rna, rnaseq)

cor_df_ref <- quant_ref_df %>%
    group_by(gene_name) %>%
    summarise(y = max(qPCR),
              r = round(cor(qPCR, rnaseq), 2),
              rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

others_1 <- ggplot(quant_ref_df, aes(rnaseq, qPCR)) +
    geom_point(size = .75) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_text(data = cor_df_ref, aes(x = -2.5, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3) +
    facet_wrap(~gene_name, scales = "free_y") +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "RNA-seq (Salmon ref)")

# hlamapper
hlamapper <- read_rds("./plot_data/hlamapper.rds")

quant_hlamapper_df <- hlamapper %>%
    left_join(nci_expression, by = c("sampleid", "gene_name")) %>%
    select(sampleid, gene_name, qPCR = m_rna, rnaseq = gene_counts)

cor_df_hlamapper <- quant_hlamapper_df %>%
    group_by(gene_name) %>%
    summarise(x = min(rnaseq),
              y = max(qPCR),
              r = round(cor(qPCR, rnaseq), 2),
              rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

others_2 <- ggplot(quant_hlamapper_df, aes(rnaseq, qPCR)) +
    geom_point(size = .75) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_text(data = cor_df_hlamapper, aes(x = -2.5, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3) +
    scale_x_continuous(labels = scientific,
                       breaks = pretty_breaks(3)) +
    facet_wrap(~gene_name, scales = "free") +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "RNA-seq (hlamapper::rna + featureCounts)")

# hla-mapper exon
quant_hlamapper_exon_df <- hlamapper %>%
    left_join(nci_expression, by = c("sampleid", "gene_name")) %>%
    select(sampleid, gene_name, qPCR = m_rna, rnaseq = b2m_norm_exon_counts)

cor_df_hlamapper_exon <- quant_hlamapper_exon_df %>%
    group_by(gene_name) %>%
    summarise(x = min(rnaseq),
              y = max(qPCR),
              r = round(cor(qPCR, rnaseq), 2),
              rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

others_3 <- ggplot(quant_hlamapper_exon_df, aes(rnaseq, qPCR)) +
    geom_point(size = .75) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_text(data = cor_df_hlamapper_exon, aes(x, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3) +
    facet_wrap(~gene_name, scales = "free") +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "RNA-seq exon estimates (hlamapper::rna + featureCounts) - B2M norm")

plot_grid(others_1, NULL, others_2, NULL, others_3, ncol = 1, 
          rel_heights = c(1, .05, 1, .05, 1),
          labels = c("A", "", "B", "", "C"))

ggsave("./plots/real_data_others.jpeg", height = 6)



## Timepoints
samples_timepoints <- read_tsv("./samples.txt", col_names = FALSE) %>%
    select(sampleid = X2, timepoint = X1) %>%
    group_by(sampleid) %>%
    filter(all(c(1, 2) %in% timepoint)) %>%
    ungroup() %>%
    arrange(sampleid, timepoint) %>%
    mutate(id = paste0(sampleid, "_t", timepoint))

quant_timepoints <- salmon_pers %>%
    filter(sampleid %in% samples_timepoints$id) %>%
    pivot_wider(names_from = timepoint, values_from = rnaseq)

cor_timepoints <- quant_timepoints %>%
    group_by(gene_name) %>%
    summarise(r = round(cor(t2, t1), 2),
              rho = round(cor(t2, t1, method = "spearman"), 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

p2 <- ggplot(quant_timepoints, aes(t1, t2)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    geom_text(data = cor_timepoints, aes(x = -2.5, y = 0.5, label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3.5) +
    facet_wrap(~gene_name) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())

p1 <- quant_timepoints %>%
    pivot_longer(t1:t2, names_to = "timepoint") %>%
    ggplot(aes(timepoint, value, color = timepoint)) +
    geom_boxplot(color = "grey35", fill = NA, outlier.color = NA) +
    geom_quasirandom(method = "smiley", show.legend = FALSE) +
    facet_wrap(~gene_name) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(x = NULL, y = "Std expression")

# Timepoints sequencing composition
salmon_pers_raw <- read_tsv("./salmon-pers/quants.bed.gz") %>%
    select(gene_name = id, samples_timepoints$id) %>%
    pivot_longer(-gene_name, names_to = "sampleid", values_to = "tpm") %>%
    select(sampleid, gene_name, tpm) %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_")

tpm_proportions_timepoints <- salmon_pers_raw %>%
    arrange(sampleid, timepoint, -tpm) %>%
    group_by(sampleid, timepoint) %>%
    mutate(prop = tpm/sum(tpm)) %>%
    ungroup() %>%
    mutate(labl = ifelse(prop >= 0.05, gene_name, "Others")) %>%
    group_by(sampleid, timepoint, gene = labl) %>%
    summarise(prop = sum(prop)) %>%
    ungroup()

gene_levels_timepoints <- tpm_proportions_timepoints %>%
    group_by(gene) %>%
    summarise(prop = mean(prop)) %>%
    ungroup() %>%
    arrange(-prop) %>%
    pull(gene)

others_prop <- tpm_proportions_timepoints %>%
    filter(timepoint == "t1", gene == "Others") %>%
    arrange(prop)

tpm_proportions_timepoints <- tpm_proportions_timepoints %>%
    mutate(sampleid = factor(sampleid, levels = others_prop$sampleid),
           gene = factor(gene, levels = gene_levels_timepoints))

mycols_timepoints <- c("grey", "#046C9A", "#D69C4E", "#ABDDDE",
                       "tomato3", "darkgreen", "grey35", "#ECCBAE")

plot_tp <- ggplot(tpm_proportions_timepoints, 
                  aes(sampleid, prop, fill = gene)) +
    geom_bar(stat = "identity", position = "fill", width = 0.95, alpha = .75,
             color = "grey85") +
    geom_text(aes(label = sub("%", "", scales::percent(prop, accuracy = 1))),
              size = 2.5, color = "grey25", fontface = "bold",
              position = position_stack(vjust = .85)) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = mycols_timepoints,
                      guide = guide_legend(direction = "horizontal")) +
    facet_wrap(~timepoint) +
    coord_flip() +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = "top") +
    guides(fill = guide_legend(override.aes = list(color = "black"))) +
    labs(y = "% of whole gene expression", fill = "Gene")

plot_grid(plot_tp, NULL, p1, p2, ncol = 1, 
          rel_heights = c(1, .05, .3, .3), labels = c("A", "", "B", "C"))

ggsave("./plots/timepoints.jpeg", height = 10, width = 8)

# Surface expression
surface_df <- quant_df %>%
    filter(gene_name == "HLA-C", !is.na(surface)) %>%
    pivot_longer(qPCR:rnaseq, names_to = "method", values_to = "mRNA")

cor_surface <- surface_df %>%
    group_by(method) %>%
    summarise(x = min(mRNA),
              y = max(surface),
              r = round(cor(surface, mRNA), 2),
              rho = round(cor(surface, mRNA, method = "spearman"), 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

ggplot(surface_df, aes(mRNA, surface)) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_point(size = .75) +
    geom_text(data = cor_surface, aes(x, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3) +
    facet_wrap(~method, scales = "free_x") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())

ggsave("./plots/surface.jpeg", width = 4, height = 2)

# HLA allele-level expression
salmon_allele <- read_rds("./plot_data/salmon_pers_allele.rds") %>%
    filter(grepl("_t1$", sampleid)) %>% 
    mutate(sampleid = sub("_t1$", "", sampleid)) %>%
    group_by(sampleid, gene_name) %>%
    mutate(zyg = ifelse(n_distinct(allele) == 1L, "hom", "het")) %>%
    ungroup() %>%
    mutate(lineage = sub("^([^:]+).*$", "\\1", allele)) %>%
    group_by(lineage) %>%
    filter(n() > 5L) %>%
    ungroup() %>%
    select(sampleid, gene_name, lineage, zyg, rna = tpm)

nci_allele <- read_rds("./plot_data/nci_allele.rds") %>%
    filter(sampleid %in% unique(salmon_allele$sampleid)) %>%
    group_by(sampleid, gene_name) %>%
    mutate(zyg = ifelse(n_distinct(allele) == 1L, "hom", "het")) %>%
    ungroup() %>%
    select(sampleid, gene_name, lineage = allele, zyg, rna) %>%
    group_by(lineage) %>%
    filter(n() > 5L) %>%
    ungroup()

allele_df <- 
    bind_rows("qPCR" = nci_allele, "RNA-seq" = salmon_allele, .id = "method")

allele_pcr <- allele_df %>%
    filter(method == "qPCR") %>%
    ggplot(aes(x = reorder(lineage, rna), y = rna)) +
    geom_quasirandom(alpha = 0.5, color = "midnightblue", method = "smiley") +
    geom_boxplot(fill = NA, outlier.color = NA) +
    facet_wrap(~gene_name, ncol = 1, scales = "free") +
    labs(x = NULL, y = NULL, title = "qPCR") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          title = element_text(size = 8))

allele_rnaseq <- allele_df %>%
    filter(method == "RNA-seq") %>%
    group_by(gene_name) %>%
    mutate(rna = scale(rna)) %>%
    ungroup() %>%
    ggplot(aes(x = reorder(lineage, rna), y = rna)) +
    geom_quasirandom(alpha = .5, color = "midnightblue", method = "smiley") +
    geom_boxplot(fill = NA, outlier.color = NA) +
    facet_wrap(~gene_name, ncol = 1, scales = "free") +
    labs(x = NULL, y = NULL, title = "RNA-seq") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          title = element_text(size = 8))

plot_grid(allele_pcr, allele_rnaseq, nrow = 1)

ggsave("./plots/allele_level_expression.jpg", width = 6, height = 5)

# investigate alleles contributing to differences

allele_pcr_v2 <- allele_df %>%
    filter(method == "qPCR", gene_name == "HLA-A") %>%
    ggplot(aes(x = reorder(lineage, rna), y = rna)) +
    geom_quasirandom(aes(color = sampleid == "66K00634"), method = "smiley") +
    scale_color_manual(values = c("TRUE" = "tomato3", "FALSE" = "grey70")) +
    geom_boxplot(fill = NA, outlier.color = NA) +
    facet_wrap(~gene_name, ncol = 1, scales = "free") +
    labs(x = NULL, y = NULL, title = "qPCR") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          title = element_text(size = 8))

allele_rnaseq_v2 <- allele_df %>%
    filter(method == "RNA-seq", gene_name == "HLA-A") %>%
    group_by(gene_name) %>%
    mutate(rna = scale(rna)) %>%
    ungroup() %>%
    ggplot(aes(x = reorder(lineage, rna), y = rna)) +
    geom_quasirandom(aes(color = sampleid == "66K00634"), method = "smiley") +
    scale_color_manual(values = c("TRUE" = "tomato3", "FALSE" = "grey70")) +
    geom_boxplot(fill = NA, outlier.color = NA) +
    facet_wrap(~gene_name, ncol = 1, scales = "free") +
    labs(x = NULL, y = NULL, color = "individual\n66K00634", title = "RNA-seq") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          title = element_text(size = 8))

leg_tmp <- get_legend(allele_rnaseq_v2)

plot_grid(allele_pcr_v2 + theme(legend.position = "none"), 
          allele_rnaseq_v2 + theme(legend.position = "none"),
          leg_tmp, 
          rel_widths = c(1, 1, .33),
          nrow = 1)

ggsave("./plots/allele_hla_a.jpg", width = 6, height = 2)
    


# rank_df <- allele_df %>%
#     group_by(method, gene_name, lineage) %>%
#     summarise(median = median(rna)) %>%
#     mutate(rank = rank(median)) %>%
#     ungroup() %>%
#     arrange(method, gene_name, rank) %>%
#     select(-median) %>%
#     group_by(lineage) %>% 
#     mutate(diff = abs(diff(rank))) %>%
#     ungroup()
# 
# 
# ggplot(rank_df, aes(method, rank, group = lineage)) +
#     geom_line(aes(color = diff)) +
#     scale_color_viridis_c(option = "inferno") +
#     geom_label(aes(label = lineage)) +
#     facet_wrap(~gene_name, ncol = 1) +
#     theme_minimal() +
#     theme(panel.grid = element_blank(),
#           axis.ticks.y = element_blank(),
#           axis.text.y = element_blank(),
#           legend.position = "bottom") +
#     labs(x = NULL, y = NULL, color = "rank difference")
# 
# 
# rank_df2 <- rank_df %>%
#     select(-diff) %>%
#     pivot_wider(names_from = method, values_from = rank) %>%
#     mutate(rank_diff = abs(qPCR - `RNA-seq`)) %>%
#     mutate(lineage = fct_reorder(lineage, rank_diff)) %>%
#     select(gene_name, lineage, rank_diff)
# 
# rank_mean_diff <- rank_df2 %>%
#     group_by(gene_name) %>%
#     summarise(m = mean(rank_diff))
# 
# ggplot(rank_df2, aes(x = lineage, y = rank_diff, fill = rank_diff)) +
#     geom_bar(stat = "identity", position = "stack") +
#     scale_fill_viridis_c(option = "cividis") +
#     facet_wrap(~gene_name, ncol = 1, scales="free_x") +
#     theme_bw() +
#     theme(axis.title.x = element_blank(),
#           panel.grid = element_blank())

# B2M

b2m <- salmon_pers_tpm %>%
    filter(id == "B2M") %>%
    select(gene_name = id, starts_with("66K")) %>%
    pivot_longer(-gene_name, names_to = "sampleid", values_to = "rna") %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_") %>%
    filter(timepoint == "t1") %>%
    select(sampleid, gene_name, rna) %>%
    pivot_wider(names_from = gene_name, values_from = rna)

b2m_hla_df <- salmon_pers_hla %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    left_join(b2m, by = "sampleid") %>%
    mutate(norm_exp = rna/B2M) %>%
    select(sampleid, gene_name, norm_exp) %>%
    left_join(nci_expression, by = c("sampleid", "gene_name"))

ggplot(b2m_hla_df, aes(norm_exp, m_rna)) +
    geom_point() +
    facet_wrap(~gene_name, scales = "free")

b2m_hla_df %>%
    arrange(desc(norm_exp))

b2m_plot1 <- b2m %>%
    pivot_longer(B2M, names_to = "gene", values_to = "TPM") %>%
    ggplot(aes(x = gene, y = TPM)) +
    geom_quasirandom(size = .8, method = "smiley") +
    scale_y_continuous(limits = c(0, NA)) +
    theme_bw() +
    labs(x = NULL)

quant_resid <- quant_df %>%
    select(-surface) %>%
    drop_na() %>%
    split(.$gene_name) %>%
    map_df(~mutate(., resid = lm(qPCR ~ rnaseq, data = .)$residuals)) %>%
    left_join(b2m)
    
quant_resid %>%
    filter(gene_name == "HLA-A") %>%
    arrange(desc(resid))

ggplot(quant_resid, aes(B2M, resid)) +
    geom_point() +
    facet_wrap(~gene_name)

b2m_df <- b2m %>%
    mutate(B2M = scale(B2M)[, 1]) %>%
    left_join(quant_df, .) 

b2m_plot2 <- ggplot(b2m_df, aes(rnaseq, qPCR)) +
    geom_point(aes(color = B2M), size = .8) +
    scale_color_gradient2(guide = guide_colourbar(barwidth = .25)) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~gene_name, scales = "free_y") +
    theme(panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey45"),
          panel.grid = element_blank()) +
    labs(x = "RNA-seq")


plot_grid(plot_grid(NULL, b2m_plot1, NULL, nrow = 1, rel_widths = c(.06, .5, 1)),
          NULL,
          plot_grid(NULL, b2m_plot2, rel_widths = c(.06, 1)),
          ncol = 1, rel_heights = c(.7, 0.1, 1), labels = c("A)", NULL, "B)"))

ggsave("./plots/b2m.jpeg", width = 6, height = 3.5)

# Concentration

salmon_pers_raw_all <- read_tsv("./salmon-pers/quants.bed.gz") %>%
    select(gene_name = id, ends_with("t1")) %>%
    pivot_longer(-gene_name, names_to = "sampleid", values_to = "tpm") %>%
    select(sampleid, gene_name, tpm) %>%
    mutate(sampleid = sub("_t1", "", sampleid))

tpm_proportions_pers_all <- salmon_pers_raw_all %>%
    arrange(sampleid, -tpm) %>%
    group_by(sampleid) %>%
    mutate(prop = tpm/sum(tpm)) %>%
    ungroup() %>%
    mutate(labl = ifelse(prop >= 0.05, gene_name, "Others")) %>%
    group_by(sampleid, gene = labl) %>%
    summarise(prop = sum(prop)) %>%
    ungroup()

concentration_df <- salmon_pers_raw_all %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    left_join(filter(tpm_proportions_pers_all, gene == "Others") %>%
                  select(sampleid, others = prop))

cor_concent <- concentration_df %>%
    group_by(gene_name) %>%
    summarise(x = min(tpm),
              y = max(others),
              r = round(cor(others, tpm), 2),
              rho = round(cor(others, tpm, method = "spearman"), 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()


ggplot(concentration_df, aes(tpm, others)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    geom_text(data = cor_concent, aes(x, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3) +
    facet_wrap(~gene_name, scales = "free_x") +
    labs(x = "TPM", y = "1 - concentration")

ggsave("./plots/hla_vs_concentration.jpeg", width = 6, height = 2)


# salmon HLA counts normalized by B2M
# salmon_b2m_norm <- read_rds("./plot_data/salmon_pers_b2m.rds") %>%
#     left_join(nci_expression) %>%
#     select(sampleid, gene_name, rnaseq = norm_counts, qPCR = m_rna) %>%
#     filter(rnaseq < 8)
# 
# cor_b2m_norm <- salmon_b2m_norm %>%
#     group_by(gene_name) %>%
#     summarise(y = max(qPCR),
#               r = round(cor(qPCR, rnaseq), 2),
#               rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
#               rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
#     ungroup()
# 
# ggplot(salmon_b2m_norm, aes(rnaseq, qPCR)) +
#     geom_point() +
#     geom_smooth(method = "lm", se = FALSE) +
#     geom_text(data = cor_b2m_norm, aes(x = -.5, y, label = rho_lab),
#               hjust = "inward", vjust = "inward", 
#               parse = TRUE,
#               size = 3.5) +
#     facet_wrap(~gene_name, scales = "free") +
#     labs(x = "RNA-seq")




    




# # remove 1 and 10th decile of B2M expression
# b2m_filtered <- b2m_df %>%
#     mutate(q = ntile(B2M, 10)) %>%
#     filter(! q %in% c(1, 10))
# 
# cor_b2m <- b2m_filtered %>%
#     group_by(gene_name) %>%
#     summarise(y = max(qPCR),
#               r = round(cor(qPCR, rnaseq), 2),
#               rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
#               rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
#     ungroup()
# 
# b2m_filtered %>%
#     filter(! sampleid %in% c("66K00046", "66K00634"))
#     
# ggplot(b2m_filtered, aes(rnaseq, qPCR)) +
#     geom_point(aes(color = B2M)) +
#     scale_color_gradient2() +
#     geom_smooth(method = "lm", se = FALSE) +
#     geom_text(data = cor_b2m, aes(x = -2, y, label = rho_lab),
#               hjust = "inward", vjust = "inward", 
#               parse = TRUE,
#               size = 3.5) +
#     facet_wrap(~gene_name, scales = "free_y") +
#     theme(panel.grid.minor = element_blank(),
#           panel.background = element_rect(fill = "grey45"),
#           panel.grid = element_blank(),
#           legend.position = "top") +
#     labs(x = "RNA-seq")
# 
# ## remove outlier individuals for HLA-A
# b2m_filtered_2 <- b2m_df %>%
#     filter(! sampleid %in% c("66K00046", "66K00634"))
# 
# cor_b2m_2 <- b2m_filtered_2 %>%
#     group_by(gene_name) %>%
#     summarise(y = max(qPCR),
#               r = round(cor(qPCR, rnaseq), 2),
#               rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
#               rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
#     ungroup()
# 
# ggplot(b2m_filtered_2, aes(rnaseq, qPCR)) +
#     geom_point(aes(color = B2M)) +
#     scale_color_gradient2() +
#     geom_smooth(method = "lm", se = FALSE) +
#     geom_text(data = cor_b2m_2, aes(x = -2, y, label = rho_lab),
#               hjust = "inward", vjust = "inward", 
#               parse = TRUE,
#               size = 3.5) +
#     facet_wrap(~gene_name, scales = "free_y") +
#     theme(panel.grid.minor = element_blank(),
#           panel.background = element_rect(fill = "grey45"),
#           panel.grid = element_blank(),
#           legend.position = "top") +
#     labs(x = "RNA-seq")
#     