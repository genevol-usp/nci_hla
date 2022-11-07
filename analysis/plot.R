library(tidyverse)
library(scales)
library(ggthemes)
library(ggbeeswarm)
library(cowplot)
library(tidytext)
library(ggrepel)
library(ggsci)
library(pals)

# Main figs

#### locus overall expression
salmon_pers_tpm <- read_tsv("./salmon-pers/quants.bed.gz") 
salmon_ref_tpm <- read_tsv("./salmon/quants.bed.gz")

salmon_pers_hla <- salmon_pers_tpm %>%
    filter(id %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    pivot_longer(starts_with("66K"), names_to = "sampleid", values_to = "rna") %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_") %>%
    filter(timepoint == "t1") %>%
    select(sampleid, gene_name = id, rna)

salmon_ref_hla <- salmon_ref_tpm %>%
    filter(id %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    pivot_longer(starts_with("66K"), names_to = "sampleid", values_to = "rna") %>%
    select(sampleid, gene_name = id, rna)

fig2_a <- salmon_pers_hla %>%
    ggplot(aes(fct_reorder(gene_name, rna, .fun = "median", .desc = TRUE), rna)) +
    geom_quasirandom(size = .75, alpha = .5, method = "smiley") +
    theme(text = element_text(size = 8, family = "Helvetica"),
          plot.margin = margin(.5, .5, .5, .5, unit = "cm"),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey96")) +
    labs(x = NULL, y = "TPM", 
         title = "RNAseq-derived HLA gene\nexpression in a cohort of\n96 individuals.")


salmon_ref_hla %>%
    mutate(gene_name = sub("HLA-", "", gene_name)) %>%
    ggplot(aes(fct_reorder(gene_name, rna, .fun = "median", .desc = TRUE), rna)) +
    geom_quasirandom(size = .25, alpha = .5, method = "smiley") +
    theme_bw() +
    theme(axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 6),
          plot.title = element_text(size = 6),
          plot.margin = margin(.5, .5, .5, .5, unit = "cm")) +
    labs(x = NULL, y = "TPM", 
         title = "RNAseq-derived HLA gene expression\n for a cohort of 96 individuals.")


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

salmon_pers_std <- salmon_pers_hla %>%
    group_by(gene_name) %>%
    mutate(rna = GenABEL::rntransform(rna)) %>%
    ungroup()

salmon_ref_std <- salmon_ref_hla %>%
    group_by(gene_name) %>%
    mutate(rna = GenABEL::rntransform(rna)) %>%
    ungroup()

# salmon_mapper_std <- read_tsv("./pipeline_results/salmon/quants.bed") %>%
#     filter(id %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
#     select(gene_name = id, starts_with("66K")) %>%
#     pivot_longer(-gene_name, names_to = "sampleid", values_to = "rnaseq") %>%
#     group_by(gene_name) %>%
#     mutate(rnaseq = GenABEL::rntransform(rnaseq)) %>%
#     ungroup()

quant_df <- 
    bind_rows("Personalized" = salmon_pers_std, 
	      "Ref transcriptome" = salmon_ref_std,
	     # "Ref hla-mapper" = salmon_mapper_std,
	      .id = "method") %>%
    left_join(nci_expression, by = c("sampleid", "gene_name")) %>%
    mutate(lab = paste0(gene_name, " (", method, ")")) %>%
    select(sampleid, gene_name, method, lab, qPCR = m_rna, rnaseq = rna, surface)

cor_df <- quant_df %>%
  group_by(lab) %>%
  summarise(y = max(qPCR),
            r = round(cor(qPCR, rnaseq), 2),
            rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
            p = cor.test(rnaseq, qPCR, method = "pearson", exact = FALSE)$p.value,
            p = format(p, format = "e", digits = 2),
            rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
  ungroup()

fig2_b <- quant_df %>%
  filter(method == "Personalized") %>%
  ggplot(aes(rnaseq, qPCR)) +
    geom_smooth(method = "lm", se = FALSE, color = "black", size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "white", size = 1) +
    geom_point(size = 1) +
    geom_text(data = cor_df %>% 
              filter(grepl("Person", lab)) %>% 
              mutate(gene_name = sub("^(\\S+).*$", "\\1", lab)), 
            aes(x = -2.5, y + (y *.1), label = rho_lab),
            hjust = "inward", vjust = "inward", 
            parse = TRUE,
            size = 2.5) +
  facet_wrap(~gene_name, nrow = 1, scales = "free") +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 8, family = "Helvetica"),
        plot.margin = margin(.5, .5, 0, .5, unit = "cm"),
        panel.background = element_rect(fill = "grey96")) +
  labs(x = "RNA-seq", title = "Correlation between qPCR and RNAseq-derived expression levels for Class I HLA genes.")

ggsave("./plots/fig2.jpg", fig2_b, width = 6, height = 2)



# Surface expression
surface_df <- quant_df %>%
    filter(method == "Personalized", 
           gene_name == "HLA-C", 
           !is.na(surface)) %>%
    select(-method) %>%
    rename("RNAseq" = "rnaseq") %>%
    pivot_longer(qPCR:RNAseq, names_to = "method", values_to = "mRNA")

cor_surface <- surface_df %>%
    group_by(method) %>%
    summarise(x = min(mRNA),
              y = max(surface),
              r = round(cor(surface, mRNA), 2),
              rho = round(cor(surface, mRNA, method = "spearman"), 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

fig2_c <- ggplot(surface_df, aes(mRNA, surface)) +
    geom_smooth(method = "lm", se = FALSE, color = "black", size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "white", size = 1) + 
    geom_point(size = 1) +
    geom_text(data = cor_surface, aes(x, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 2.5) +
    facet_wrap(~method, scales = "free_x", strip.position = "bottom") +
    theme(strip.background = element_blank(),
          strip.placement = "outside") +
    theme(panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey96"),
          text = element_text(size = 8, family = "Helvetica"),
          strip.text = element_text(vjust = 3),
          plot.margin = margin(.5, .5, 0, .5, unit = "cm")) +
    labs(x = NULL, y = "Surface levels",
         title = "Correlation between HLA-C surface protein\nlevels with qPCR or RNAseq-derived mRNA\nlevels in 48 individuals.")


fig2 <- plot_grid(
    plot_grid(fig2_a, fig2_c, nrow = 1, rel_widths = c(.7, 1),
              labels = c("A)", "C)"), label_size = 10, label_fontfamily = "Helvetica"),
    fig2_b, labels = c(NA, "B)"), label_size = 10,
    ncol = 1)

ggsave("./plots/fig2.jpg", fig2, width = 6, height = 4)



# HLA allele-level expression
salmon_allele <- read_rds("./plot_data/salmon_pers_allele.rds") %>%
    filter(grepl("_t1$", sampleid)) %>% 
    mutate(sampleid = sub("_t1$", "", sampleid),
           lineage = sub("^([^:]+).*$", "\\1", allele)) %>%
    select(sampleid, gene_name, lineage, rna = tpm)

nci_allele <- read_rds("./plot_data/nci_allele.rds") %>%
    filter(sampleid %in% unique(salmon_allele$sampleid)) %>%
    select(sampleid, gene_name, lineage = allele, rna, rna_adj)

alleles_filt <-  
    bind_rows("qPCR" = nci_allele, "RNA-seq" = salmon_allele, .id = "method") %>%
    split(.$method) %>%
    map(~group_by(., lineage) %>% 
            filter(n() >= 5) %>%
            ungroup()) %>%
    map("lineage") %>%
    reduce(.f = intersect) %>%
    sort()

allele_df <- 
    bind_rows("qPCR" = nci_allele, "RNA-seq" = salmon_allele, .id = "method") %>%
    filter(lineage %in% alleles_filt)

plot_ranks <- function(gene, allele_data) {
    
    plot_df_a <- allele_data %>%
        filter(gene_name == gene) %>%
        group_by(method, gene_name, lineage) %>%
        summarise(m = median(rna)) %>%
        ungroup() %>%
        split(.$method) %>%
        map_df(~arrange(., m) %>% mutate(index = row_number())) %>%
        select(-m) 
    
    a_diffs <- plot_df_a %>%
        pivot_wider(names_from = method, values_from = index) %>%
        mutate(d = abs(qPCR - `RNA-seq`)) %>%
        select(gene_name, lineage, d)
    
    rank_a <- plot_df_a %>%
        left_join(a_diffs) %>%
        mutate(y = as.numeric(factor(method, levels = c("qPCR", "RNA-seq"))),
               y2 = case_when(y == 1L ~ y + .15,
                              y == 2L ~ y - .15)) %>%
        ggplot(aes(x = index, y = y, label = lineage)) +
        geom_line(aes(x = index, y = y2, group = lineage), size = .2, color = "grey60") +
        geom_label(aes(y = y, fill = d), 
                   size = 2.5, 
                   fontface = "bold",
                   label.padding = unit(0.1, "lines"),
                   show.legend = FALSE) +
        scale_fill_gradient(low = "white", high = "grey50") +
        scale_y_continuous(breaks = c(1L, 2L), 
                           limits = c(.8, 2.2),
                           labels = c("qPCR", "RNA-seq")) +
        theme_minimal() +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 8, family = "Helvetica"),
              axis.title = element_blank(),
              panel.grid = element_blank())
    
    qpcr_a <- allele_data %>%
        filter(gene_name == gene, method == "qPCR") %>%
        ggplot(aes(x = reorder(lineage, rna, "median"), y = rna)) +
        geom_boxplot(fill = NA, alpha = .25, size = .25, color = "grey50",
                     outlier.color = NA) +
        geom_quasirandom(method = "smiley", fill = "white", shape = 21,
                         alpha = .75, size = 1.5, width = .25) +
        scale_y_reverse(labels = function(x) format(x, nsmall = 1)) +
        theme_minimal() +
        theme(axis.title = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 8, family = "Helvetica"),
              panel.grid = element_blank(),
              plot.margin = margin(l = 25, t = -10, b = 10))
    
    rnaseq_a <- allele_data %>%
        filter(gene_name == gene, method == "RNA-seq") %>%
        ggplot(aes(x = reorder(lineage, rna, "median"), y = rna/1000L)) +
        geom_boxplot(fill = NA, alpha = .25, size = .25, color = "grey50",
                     outlier.color = NA) +
        geom_quasirandom(method = "smiley", fill = "black", shape = 21,
                         alpha = .75, size = 1.5, width = .25) +
        scale_y_continuous(labels = function(x) format(x, nsmall = 1)) +
        theme_minimal() +
        theme(axis.title = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 8, family = "Helvetica"),
              panel.grid = element_blank(),
              plot.margin = margin(l = 25, b = -10, t = 10))
    
    plot_grid(rnaseq_a, rank_a, qpcr_a, ncol = 1, rel_heights = c(1, 1.25, 1))
}

fig3 <- plot_grid(plot_ranks("HLA-A", allele_df), 
                  plot_ranks("HLA-B", allele_df), 
                  plot_ranks("HLA-C", allele_df),
                  ncol = 1, 
                  labels = c("HLA-A", "HLA-B", "HLA-C"),
                  label_fontfamily = "Helvetica", label_size = 10,
                  hjust = 0, label_x = .5) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/fig3.jpg", fig3, height = 6.25, width = 5)

ggsave("./plots/fig3_adj.jpg", fig3, height = 6.25, width = 5)



# Fig 4
# 4A
ids <- read_lines("./sample_ids.txt")

tx_annots <- read_tsv("../indices/transcript_annotation_df.tsv") %>%
    select(tx_id, gene_id, gene_name, tx_type)

salmon_tx <- "./salmon/quant/%s/quant.sf" %>%
    sprintf(ids) %>%
    setNames(ids) %>%
    map_df(. %>% read_tsv %>% select(tx_id = Name, len = EffectiveLength, counts = NumReads, tpm = TPM), 
           .id = "sampleid") %>%
    mutate(timepoint = ifelse(grepl("t1$", sampleid), "Timepoint 1", "Timepoint 2")) %>%
    left_join(tx_annots) %>%
    select(sampleid, timepoint, tx_id, tx_type, len, gene_id, gene_name, counts, tpm)

samples_timeps <- tibble(x = read_lines("./sample_ids.txt")) %>%
    separate(x, c("sampleid", "timepoint"), sep = "_", remove = FALSE) %>%
    group_by(sampleid) %>%
    filter(all(c("t1", "t2") %in% timepoint)) %>%
    ungroup() %>%
    select(sampleid = x, x = sampleid, timepoint)

quant_batches_tpm <- salmon_tx %>%
    filter(sampleid %in% samples_timeps$sampleid) %>%
    filter(tx_type == "protein_coding") %>%
    group_by(sampleid, timepoint, gene_id, gene_name) %>%
    summarise(tpm = sum(tpm)) %>%
    ungroup() %>%
    mutate(timepoint = sub("imepoint ", "", timepoint),
           sampleid = sub("_t[12]$", "", sampleid)) %>%
    pivot_wider(names_from = timepoint, values_from = tpm)

expressed_genes <- quant_batches_tpm %>%
    group_by(gene_id, gene_name) %>%
    filter(mean(T1 > 25) > .5) %>%
    ungroup()

cor_tps_df <- expressed_genes %>%
    group_by(sampleid) %>%
    summarise(rho = round(cor(T2, T1, method = "spearman"), 2),
              rho_lab = paste("~rho == ", rho)) %>%
    ungroup()

fig4_a <- ggplot(expressed_genes, aes(log10(T1+1), log10(T2+1))) +
    geom_abline(linetype = 2) +
    geom_point(size = .25, alpha = .25) +
    geom_text(data = cor_tps_df, 
              aes(x = 0, y = 4, label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3) +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    facet_wrap(~sampleid, nrow = 2) +
    theme(text = element_text(size = 10, family = "Helvetica"),
          panel.background = element_rect(fill = "grey96"),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = "white"),
          plot.title = element_text(size = 9),
          plot.margin = margin(.5, .5, .5, .5, unit = "cm")) +
    labs(x = "Original sample", y = "Fresh sample") +
    labs(x = "Original sample", y = "Fresh sample",
         title = "High correlation between two RNA-seq experiments for 11 individuals.")


alleles_df <- read_rds("./plot_data/salmon_pers_allele.rds") %>%
    filter(grepl("_t1$", sampleid)) %>% 
    mutate(sampleid = sub("_t1$", "", sampleid),
           lineage = sub("^([^:]+).*$", "\\1", allele)) %>%
    select(sampleid, gene_name, lineage) %>%
    group_by(sampleid, gene_name) %>%
    mutate(a = paste0("a", 1:n())) %>%
    ungroup() %>%
    pivot_wider(names_from = a, values_from = lineage)

# dose_df <- quant_df %>%
#     filter(method == first(method)) %>%
#     select(-method, -surface, -lab) %>%
#     left_join(alleles_df, by = c("sampleid", "gene_name")) %>%
#     pivot_longer(a1:a2, names_to = "hap", values_to = "alleles") %>%
#     group_by(alleles) %>%
#     filter(n_distinct(sampleid) >= 15) %>%
#     ungroup() %>%
#     mutate(dose = 1) %>%
#     group_by(sampleid, gene_name, qPCR, rnaseq, alleles) %>%
#     summarise(dose = sum(dose)) %>%
#     ungroup() %>%
#     mutate(alleles = factor(alleles, levels = sort(unique(alleles)))) %>%
#     filter(dose > 0)
# 
# fig4_b <- ggplot(dose_df, aes(rnaseq, qPCR)) +
#     geom_line(stat = "smooth", method = "loess", span = 1, se = FALSE, 
#               alpha = .5, size = 1, color = "black") +
#     geom_point(aes(fill = factor(dose)), size = 1.5, stroke = .25, shape = 21) +
#     scale_fill_manual(values = c("1" = "grey50", "2" = "black")) +
#     scale_x_continuous(breaks = pretty_breaks(3)) +
#     scale_y_continuous(breaks = pretty_breaks(3)) +
#     facet_wrap(~alleles, scales="free") +
#     theme(text = element_text(size = 10, family = "Helvetica"),
#           panel.grid = element_blank(),
#           panel.background = element_rect(fill = "grey96"),
#           plot.margin = margin(.5, 0, 0, .5, unit = "cm"),
#           plot.title = element_text(size = 8)) +
#     labs(x = "RNA-seq", fill = "allele\ndose",
#          title = "RNA-seq vs. qPCR relationship broken by specific alleles.") +
#     guides(fill = guide_legend(override.aes = list(size = 2)))

# 4B
alleles_quant_df <- quant_df %>%
    filter(method == first(method), gene_name %in% c("HLA-A", "HLA-C")) %>%
    select(-surface, -method, -lab) %>%
    left_join(alleles_df, by = c("sampleid", "gene_name")) %>%
    mutate(`A*03` = as.integer(a1 == "A*03") + as.integer(a2 == "A*03"),
           `C*07` = as.integer(a1 == "C*07") + as.integer(a2 == "C*07"),
           gene_name = case_when(gene_name == "HLA-A" ~ "A*03",
                                 gene_name == "HLA-C" ~ "C*07",
                                 TRUE ~ NA_character_),
           dose = case_when(gene_name == "A*03" ~ `A*03`,
                            gene_name == "C*07" ~ `C*07`)) %>%
    select(sampleid, allele = gene_name, qPCR, rnaseq, dose)

fig4_b <- ggplot(alleles_quant_df, aes(rnaseq, qPCR)) +
    geom_point(aes(fill = factor(dose)), 
               shape = 21, stroke = .2, size = 1.5, color = "black") +
    scale_fill_manual(values = c("0" = "white", "1" = "grey50", "2" = "black")) +
    facet_wrap(~allele, scale = "free") +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size = 10, family = "Helvetica"),
          panel.background = element_rect(fill = "grey96"),
          legend.key.width = unit(0.25, "cm"),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-5,-5,-5,-5),
          plot.margin = margin(.25, .5, 0, .5, unit = "cm"),
          plot.title = element_text(size = 9)) +
    labs(x = "RNA-seq", fill = "allele\ndose",
         title = "Correlation RNA-seq vs qPCR highlighting\nthe number of copies of A*03 and C*07")
    

#4C
# new a01 a11 expression
newadf <-
    read_tsv("./plot_data/a01a11_quants.tsv") %>%
    group_by(sampleid, gene_name) %>%
    summarise(rna = sum(tpm)) %>%
    ungroup()

newa_std <- anti_join(salmon_pers_hla, distinct(newadf, sampleid, gene_name)) %>%
    bind_rows(newadf) %>%
    filter(gene_name == "HLA-A") %>%
    group_by(gene_name) %>%
    mutate(rna = GenABEL::rntransform(rna)) %>%
    ungroup()

quant_newa_df <- 
    left_join(newa_std, nci_expression, by = c("sampleid", "gene_name")) %>%
    select(sampleid, gene_name, qPCR = m_rna, rnaseq = rna)

cor_newa_df <- quant_newa_df %>%
    group_by(gene_name) %>%
    summarise(y = max(qPCR),
              r = round(cor(qPCR, rnaseq), 2),
              rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
              p = cor.test(rnaseq, qPCR, method = "pearson", exact = FALSE)$p.value,
              p = format(p, format = "e", digits = 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

fig4_c <- quant_newa_df %>%
    ggplot(aes(rnaseq, qPCR)) +
    geom_smooth(method = "lm", se = FALSE, color = "black", size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "white", size = 1) +
    geom_point(size = 1) +
    geom_text(data = cor_newa_df,
              aes(x = -2.5, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3) +
    facet_wrap(~gene_name, nrow = 1, scales = "free") +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size = 10, family = "Helvetica"),
          plot.margin = margin(.25, .5, 0, 0, unit = "cm"),
          plot.title = element_text(size = 9),
          panel.background = element_rect(fill = "grey96")) +
    labs(x = "RNA-seq", title = "Correlation after ajusting for\nisoform length for A*01 and A*11")


#4D
samples_t1 <- grep("_t1$", ids, value = TRUE)

salmon_pers_tx <- "./salmon-pers/quant/%s/quant.sf" %>%
    sprintf(samples_t1) %>%
    setNames(samples_t1) %>%
    map_df(. %>% read_tsv %>% select(tx_id = Name, counts = NumReads, tpm = TPM), 
           .id = "sampleid") %>%
    mutate(tx_id = sub("^([^_]+).*$", "\\1", tx_id),
	   sampleid = sub("_t1$", "", sampleid)) %>%
    left_join(tx_annots) %>%
    select(sampleid, tx_id, tx_type, gene_id, gene_name, tpm)

b2m_all <- salmon_pers_tx %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C", "B2M")) %>%
    group_by(sampleid, gene_name) %>%
    summarise(tpm = sum(tpm)) %>%
    ungroup() %>%
    pivot_wider(names_from = gene_name, values_from = tpm) %>%
    mutate_at(vars(-sampleid), GenABEL::rntransform) %>%
    pivot_longer(starts_with("HLA"), names_to = "gene_name", values_to = "std")

b2m_all_cor <- b2m_all %>%
    group_by(gene_name) %>%
    summarise(x = min(std),
	      y = max(B2M),
              rho = round(cor(B2M, std, method = "spearman"), 2),
              rho_lab = paste("~rho == ", rho)) %>%
    ungroup()

fig4_d <- ggplot(b2m_all, aes(std, B2M)) +
    geom_smooth(method = "lm", se = FALSE, color = "black", size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "white", size = 1) +
    geom_point() +
    geom_text(data = b2m_all_cor, 
	      aes(x = x, y = y, label = rho_lab),
	      hjust = "inward", vjust = "inward", 
	      parse = TRUE,
	      size = 3) +
    facet_wrap(~gene_name, nrow = 1, scales = "free") +
    theme(panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey96"),
          text = element_text(size = 10, family = "Helvetica"),
          plot.margin = margin(.5, 1, .5, .5, unit = "cm"),
          plot.title = element_text(size = 9)) +
    labs(x = "HLA Std-normal expression",
         title = "Correlation between HLA and B2M expression.")


fig4 <- 
    plot_grid(fig4_a,
              plot_grid(fig4_b, fig4_c, rel_widths = c(1, .6),
                        labels = c("B)", "C)"), 
                        label_size = 10, label_fontfamily = "Helvetica"),
              fig4_d, 
              ncol = 1, 
              labels = c("A)", NA, "D)"),
              label_size = 10, label_fontfamily = "Helvetica",
              rel_heights = c(1, .7, .7))


ggsave("./plots/fig4.jpg", fig4, width = 6, height = 7)






###############################################################################
# Supplements




## Fig S2

# Read coverage
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

fig_s1_a <- ggplot(reads_t1, aes(reorder(sampleid, depth), depth)) +
    geom_col(width = 1) +
    scale_y_continuous(labels = function(x) round(floor(x)/1e6),
                       breaks = pretty_breaks(10)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.margin = margin(.5, .5, .5, .5, unit = "cm"),
          plot.title = element_text(size = 10)) +
    labs(x = NULL, y = "Million reads", 
         title = "Depth of sequencing per individual.")

# Biotypes
meantpm_by_types <- salmon_tx %>%
    group_by(timepoint) %>%
    mutate(pct_tpm = tpm/sum(tpm)) %>%
    group_by(timepoint, tx_type) %>%
    summarise(pct_tpm = sum(pct_tpm)) %>%
    mutate(tx_type_lab = ifelse(pct_tpm > 0.01, tx_type, "Other")) %>%
    group_by(timepoint, tx_type_lab) %>%
    summarise(pct_tpm = sum(pct_tpm)) %>%
    ungroup() %>%
    arrange(timepoint, desc(pct_tpm)) %>%
    mutate(tx_type_lab = recode(tx_type_lab, "nonsense_mediated_decay" = "NMD"),
           tx_type_lab = fct_inorder(tx_type_lab),
           tx_type_lab = fct_rev(tx_type_lab))

cols <- c("grey", pal_npg()(10), "black")

fig_s1_b <- meantpm_by_types %>%
    filter(timepoint == "Timepoint 1") %>%
    ggplot(aes(tx_type_lab, pct_tpm, fill = tx_type_lab)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = cols) +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.margin = margin(.5, .5, 0, .5, unit = "cm"),
          plot.title = element_text(size = 10, margin = margin(b = -10))) +
    labs(x = NULL, y = NULL,
         title = "Proportion of transcripts detected for\neach Gencode biotype.")

# PCA
quant_all_tpm <- salmon_pers_tpm %>%
    mutate(gid = sub("\\.\\d+$", "", gid)) %>%
    select(gene_id = gid, gene_name = id, starts_with("66K")) %>%
    pivot_longer(-(gene_id:gene_name), names_to = "sampleid") %>%
    filter(grepl("_t1$", sampleid)) %>%
    mutate(sampleid = sub("_t1$", "", sampleid)) %>%
    group_by(gene_id, gene_name) %>%
    filter(mean(value > 1) > 0.5) %>%
    ungroup()

# var_genes <- salmon_pers_tpm %>%
#     mutate(gid = sub("\\.\\d+$", "", gid)) %>%
#     select(gene_id = gid, gene_name = id, starts_with("66K")) %>%
#     pivot_longer(-(gene_id:gene_name), names_to = "sampleid") %>%
#     filter(grepl("_t1$", sampleid)) %>%
#     mutate(sampleid = sub("_t1$", "", sampleid)) %>%
#     group_by(gene_id, gene_name) %>%
#     summarise(v = var(value)) %>%
#     ungroup() %>%
#     top_n(2000, v)
# 
# quant_var_genes <- salmon_pers_tpm %>%
#     mutate(gid = sub("\\.\\d+$", "", gid)) %>%
#     select(gene_id = gid, gene_name = id, starts_with("66K")) %>%
#     inner_join(var_genes) %>%
#     pivot_longer(-(gene_id:gene_name), names_to = "sampleid") %>%
#     filter(grepl("_t1$", sampleid)) %>%
#     mutate(sampleid = sub("_t1$", "", sampleid))

quant_matrix <- quant_all_tpm %>%
    select(-gene_name) %>%
    pivot_wider(names_from = gene_id, values_from = value) %>%
    column_to_rownames("sampleid") %>%
    as.matrix()

gene_ids <- distinct(quant_all_tpm, gene_id, gene_name)

pca <- prcomp(quant_matrix, center = TRUE, scale. = TRUE, rank. = 100)
pc_scores <- as_tibble(pca$x, rownames = "sampleid")

pc_loadings <- as_tibble(pca$rotation, rownames = "gene_id") %>%
    pivot_longer(-gene_id, names_to = "pc") %>%
    mutate(s = ifelse(value > 0, 1, 0)) %>%
    group_by(pc, s) %>%
    slice_max(n = 25, order_by = abs(value)) %>%
    ungroup() %>%
    left_join(gene_ids, by = "gene_id") %>%
    select(gene_id, gene_name, pc, value)

pc_loadings %>%
    filter(pc %in% c("PC1", "PC2")) %>%
    ggplot(aes(x = value, y = reorder_within(gene_name, by = value, within = pc))) +
    geom_col(aes(fill = value > 0), show.legend = FALSE) +
    scale_y_discrete(labels = function(x) str_remove(x, "_+PC\\d+$")) +
    facet_wrap(~pc, scales = "free") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(size = 10)) +
    labs(x = "PC loadings", y = NULL,
         title = "Top 25 contributions to each direction")

#ggsave("./plots/pc_loadings.png", width = 6, height = 5)


ribo <- quant_all_tpm %>%
    mutate(ribo = grepl("^RPS|^RPL", gene_name)) %>%
    group_by(sampleid, ribo) %>%
    summarise(value = sum(value)) %>%
    group_by(sampleid) %>%
    mutate(value = value/sum(value)) %>%
    ungroup() %>%
    filter(ribo == TRUE) %>%
    select(sampleid, ribo_pct = value)

srp <- quant_all_tpm %>%
    mutate(srp = grepl("^RN7S", gene_name)) %>%
    group_by(sampleid, srp) %>%
    summarise(value = sum(value)) %>%
    group_by(sampleid) %>%
    mutate(value = value/sum(value)) %>%
    ungroup() %>%
    filter(srp == TRUE) %>%
    select(sampleid, srp_pct = value)

mito <- quant_all_tpm %>%
    mutate(mito = grepl("^MT-", gene_name)) %>%
    group_by(sampleid, mito) %>%
    summarise(value = sum(value)) %>%
    group_by(sampleid) %>%
    mutate(value = value/sum(value)) %>%
    ungroup() %>%
    filter(mito == TRUE) %>%
    select(sampleid, mito_pct = value)

## just check if samples cluster according to sequencing depth
pc_scores %>%
    select(sampleid, PC1:PC2) %>%
    left_join(reads_t1) %>%
    ggplot(aes(PC1, PC2)) +
    geom_point(aes(fill = log10(depth)), size = 1.5, shape = 21) +
    scale_fill_viridis_c("", option = "magma",
                         guide = guide_colorbar(barwidth = .25,
                                                barheight = 7)) + 
    theme_bw()
##

fig_s1_c <- pc_scores %>%
    select(sampleid, PC1:PC2) %>%
    left_join(ribo) %>%
    left_join(srp) %>%
    left_join(mito) %>%
    pivot_longer(ribo_pct:mito_pct) %>%
    mutate(name = recode(name, 
                         "ribo_pct" = "% Ribosome", 
                         "srp_pct" = "% SRP",
                         "mito_pct" = "% Mitochondria")) %>%
    split(.$name) %>%
    map(~ggplot(., aes(PC1, PC2)) +
            geom_point(aes(fill = value), size = 1.5, shape = 21) +
            scale_fill_viridis_c("", option = "magma",
                                 labels = scales::percent,
                                 guide = guide_colorbar(barwidth = .25,
                                                        barheight = 5)) + 
            facet_wrap(~name, nrow = 1) +
            theme_bw() +
            theme(legend.text = element_text(size = 8),
                  legend.margin = margin(0, 0, 10, -10),
                  axis.title = element_text(size = 8),
                  axis.text = element_text(size = 6))) %>%
    plot_grid(plotlist = ., nrow = 1) +
    theme(plot.title = element_text(size = 10),
          plot.margin = margin(.5, .5, .5, .5, unit = "cm")) +
    labs(title = "Percentage of transcripts from Mitochondria, Ribosomes, and Signal Recognition Particles")
    

## PCA correction
quants_pca <- "./salmon-pers/pca_correction/corrected/quants_%s.bed.gz" %>%
    sprintf(c(1:2, seq(0, 25, 5))) %>%
    setNames(c(1:2, seq(0, 25, 5))) %>%
    map_df(~read_tsv(.) %>%
               filter(id %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
               pivot_longer(starts_with("66K"), names_to = "sampleid", values_to = "rna") %>%
               select(sampleid, gene_name = id, rna), .id = "PCs")

cor_pca <- quants_pca %>%
    mutate(PCs = factor(PCs, levels = sort(c(1:2, seq(0, 25, 5))))) %>%
    left_join(nci_expression) %>%
    select(PCs, sampleid, gene_name, rnaseq = rna, qPCR = m_rna) %>%
    group_by(PCs, gene_name) %>%
    summarise(r = round(cor(qPCR, rnaseq), 2),
              rho = round(cor(qPCR, rnaseq, method = "spearman"), 2)) %>%
    ungroup()


fig_s1_d <- ggplot(cor_pca, aes(PCs, rho, group = gene_name)) +
    geom_line(size = 2, color = "grey25", alpha = .6) +
    geom_point(size = 2, color = "grey25") +
    geom_label_repel(data = filter(cor_pca, PCs == 25),
                     aes(label = gene_name)) +
    scale_y_continuous(breaks = seq(0, .6, .1)) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.margin = margin(.5, 2.5, .5, 2.5, unit = "cm"),
          plot.title = element_text(size = 10)) +
    labs(x = "# of PCs for RNA-seq correction",
         y = "Spearman correlation\nwith qPCR estimates",
         title = "PCA correction of RNAseq-derived HLA expression estimates\ndoes not improve correlation with qPCR estimates.")


fig_s1 <- plot_grid(
    plot_grid(fig_s1_a, fig_s1_b, labels = c("A)", "B)"), label_size = 10, nrow = 1),
    fig_s1_c,
    fig_s1_d,
    labels = c(NA, "C)", "D)"), label_size = 10, ncol = 1) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/fig_s1.jpg", fig_s1, height = 7, width = 6.5)

# Fig S2
s2 <- quant_df %>%
    filter(method == "Ref transcriptome") %>%
    ggplot(aes(rnaseq, qPCR)) +
    geom_smooth(method = "lm", se = FALSE, color = "black", size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "white", size = 1) +
    geom_point(size = 1) +
    geom_text(data = cor_df %>% 
                  filter(grepl("Ref transcrip", lab)) %>% 
                  mutate(gene_name = sub("^(\\S+).*$", "\\1", lab)), 
              aes(x = -2.5, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3) +
    facet_wrap(~gene_name, nrow = 1, scales = "free") +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size = 8, family = "Helvetica"),
          plot.margin = margin(.5, .5, 0, .5, unit = "cm"),
          panel.background = element_rect(fill = "grey96")) +
    labs(x = "RNA-seq", title = "Correlation between RNA-seq and qPCR\nusing the Reference transcriptome to index RNA-seq data")

ggsave("./plots/fig_s2.jpg", s2, height = 2.5, width = 6)







# Fig S3
# New version of Salmon, with and without bias correction
quant_biascorr <- file.path("./salmon-pers/quant-bias", ids, "quant.sf") %>%
    setNames(ids) %>%
    .[grepl("_t1", .)] %>%
    map_df(. %>% read_tsv() %>%
               filter(grepl("_[ABC]\\*", Name)) %>%
               separate(Name, c("tx_id", "hla"), sep = "_") %>%
               left_join(tx_annots, by = "tx_id") %>%
               group_by(gene_id, gene_name, hla) %>%
               summarise(tpm = sum(TPM)) %>%
               group_by(gene_id, gene_name) %>%
               mutate(i = ifelse(n() == 1, "1_1", "1"),
                      tpm = ifelse(i == "1_1", tpm/2L, tpm)) %>%
               ungroup() %>%
               separate_rows(i, sep = "_") %>%
               select(-i), .id = "sampleid") %>%
    mutate(sampleid = sub("_t1", "", sampleid))

quant_nocorr <- file.path("./salmon-pers/quant-nogc", ids, "quant.sf") %>%
    setNames(ids) %>%
    .[grepl("_t1", .)] %>%
    map_df(. %>% read_tsv() %>%
               filter(grepl("_[ABC]\\*", Name)) %>%
               separate(Name, c("tx_id", "hla"), sep = "_") %>%
               left_join(tx_annots, by = "tx_id") %>%
               group_by(gene_id, gene_name, hla) %>%
               summarise(tpm = sum(TPM)) %>%
               group_by(gene_id, gene_name) %>%
               mutate(i = ifelse(n() == 1, "1_1", "1"),
                      tpm = ifelse(i == "1_1", tpm/2L, tpm)) %>%
               ungroup() %>%
               separate_rows(i, sep = "_") %>%
               select(-i), .id = "sampleid") %>%
    mutate(sampleid = sub("_t1", "", sampleid))

new_alleles_df <-
    bind_rows("Bias correction" = quant_biascorr,
              "No bias correction" = quant_nocorr, 
              .id = "method") %>%
    mutate(lineage = sub("^([ABC]\\*\\d+).*$", "\\1", hla))

alleles_order <- new_alleles_df %>%
    filter(lineage %in% alleles_filt, method == "Bias correction") %>%
    group_by(gene_name, lineage) %>%
    summarise(md = median(tpm)) %>%
    ungroup() %>%
    arrange(gene_name, md) %>%
    pull(lineage)


new_alleles_df %>%
    filter(lineage %in% alleles_filt) %>%
    mutate(lineage = factor(lineage, levels = alleles_order)) %>%
    ggplot(aes(x = lineage, y = tpm, fill = method, color = method)) +
    geom_boxplot(fill = NA, alpha = .25, size = .25,
                 outlier.color = NA, position = "dodge",
                 show.legend = FALSE) +
    geom_jitter(position = position_jitterdodge(), size = .5) +
    scale_color_npg() +
    facet_wrap(~gene_name, ncol = 1, scales = "free") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()) +
    labs(x = NULL, y = "Transcripts per Million") +
    guides(color = guide_legend(override.aes = list(size = 2)))


quant_nocorr_genes <- quant_nocorr %>%
    group_by(sampleid, gene_name) %>%
    summarise(tpm = sum(tpm)) %>%
    group_by(gene_name) %>%
    mutate(rnaseq = GenABEL::rntransform(tpm)) %>%
    ungroup()

quant_biascorr_genes <- quant_biascorr %>%
    group_by(sampleid, gene_name) %>%
    summarise(tpm = sum(tpm)) %>%
    group_by(gene_name) %>%
    mutate(rnaseq = GenABEL::rntransform(tpm)) %>%
    ungroup()

quant_bias_df <- 
    bind_rows("Bias correction" = quant_biascorr_genes,
              "No correction" = quant_nocorr_genes,
              .id = "method") %>%
    inner_join(nci_expression) %>%
    select(method, sampleid, gene_name, rnaseq, qPCR = m_rna)

cor_bias_df <- quant_bias_df %>%
    group_by(method, gene_name) %>%
    summarise(y = max(qPCR),
              r = round(cor(qPCR, rnaseq), 2),
              rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
              p = cor.test(rnaseq, qPCR, method = "pearson", exact = FALSE)$p.value,
              p = format(p, format = "e", digits = 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

fig_s3 <- quant_bias_df %>%
    filter(method == "No correction") %>%
    ggplot(aes(rnaseq, qPCR)) +
    geom_smooth(method = "lm", se = FALSE, color = "black", size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "white", size = 1) +
    geom_point(size = 1) +
    geom_text(data = cor_bias_df %>% filter(method == "No correction"), 
              aes(x = -2.5, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 2.5) +
    facet_wrap(~gene_name, nrow = 1, scales = "free") +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size = 8, family = "Helvetica"),
          plot.margin = margin(.5, .5, 0, .5, unit = "cm"),
          panel.background = element_rect(fill = "grey96")) +
    labs(x = "RNA-seq", title = "Correlation between qPCR and RNAseq-derived expression levels for Class I HLA genes\n(RNA-seq data not corrected for biases in Salmon)")

ggsave("./plots/fig_s3.jpg", fig_s3, width = 6, height = 2.5)


# Fig S4
annots <- read_tsv("../indices/transcript_annotation_df.tsv")

hla_annot <- read_rds("./plot_data/transcripts_annot.rds") %>%
    left_join(distinct(annots, transcript_id = tx_id, tx_type)) %>%
    distinct(gene_name, transcript_id, tx_type, feature, start, end)

tmp <- hla_annot %>%
    filter(tx_type == "protein_coding") %>%
    mutate(pos = map2(start, end, `:`)) %>%
    select(-start, -end) %>%
    unnest(cols = pos) %>%
    distinct(gene_name, feature, tx_type, pos) %>%
    group_by(gene_name, pos) %>%
    mutate(is_exon = ifelse(any(feature == "exon"), "Yes", "No")) %>%
    ungroup() %>%
    distinct(gene_name, pos, is_exon) %>%
    arrange(gene_name, pos)


sampleids_t1 <- read_lines("./sample_ids_t1.txt")

cov_df <- 
    paste0("./plot_data/coverage_", sampleids_t1, ".txt") %>%
    setNames(sampleids_t1) %>%
    map_df(~read_tsv(., col_names = FALSE) %>%
               select(pos = X2, cov = X3),
           .id = "sampleid")

# separate A*01 and A*11
genos <- read_tsv("./genos_final.tsv") %>%
    select(sampleid, gene_name, i, allele) %>%
    filter(gene_name == "HLA-A") %>%
    mutate(lineage = sub("^(A\\*\\d+).*$", "\\1", allele)) %>%
    group_by(sampleid) %>%
    summarise(a01_a11 = sum(lineage %in% c("A*01", "A*11"))) %>%
    ungroup() %>%
    mutate(a01_a11 = case_when(a01_a11 == 0 ~ "0 copies of A*01 and/or A*11",
                               a01_a11 == 1 ~ "1 copy of A*01 and/or A*11",
                               a01_a11 == 2 ~ "2 copies of A*01 and/or A*11"),
           a01_a11 = fct_inorder(a01_a11))

cov_df %>%
    filter(pos < 31e6) %>%
    left_join(genos) %>%
    group_by(pos, a01_a11) %>%
    summarise(cov = mean(cov)) %>%
    ungroup() %>%
    pivot_wider(names_from = a01_a11, values_from = cov) %>%
    select(pos, a0 = `0`, a2 = `2`) %>%
    tail(500) %>% print(n = Inf)


fig_s4_a <- cov_df %>%
    filter(pos < 31e6) %>%
    left_join(genos) %>%
    left_join(tmp, by = "pos") %>%
    ggplot(aes(pos, cov, color = is_exon)) +
    geom_line(aes(group = sampleid)) +
    geom_vline(xintercept = 29945765, linetype = 2) +
    scale_x_continuous(labels = function(x) round(x/1e6, 3)) +
    scale_y_continuous(labels = scales::comma) +
    scale_color_manual(values = c("Yes" = "midnightblue", "No" = "tomato3")) +
    facet_wrap(~a01_a11, ncol = 1) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = "Position in chr6 (Mb)", y = "Coverage", 
         color = "Exon:")

# test new a01 a11 expression

samples01 <- read_lines("./plot_data/a01_a11_individuals.txt")

newa1 <- read_tsv("./plot_data/a01a11_quants.tsv") %>%
    group_by(sampleid, gene_name) %>%
    mutate(n = n(),
           tpm = ifelse(n == 1, tpm/2L, tpm),
           i = ifelse(n == 1, "1_1", "1")) %>%
    ungroup()  %>%
    separate_rows(i, sep = "_") %>%
    select(-n, -i) %>%
    mutate(lineage = sub("^([ABC]\\*\\d+).*$", "\\1", allele)) %>%
    mutate(method = "RNA-seq") %>%
    select(method, sampleid, gene_name, lineage, rna = tpm)

allele_newa_df <- anti_join(allele_df, 
                       distinct(newa1, method, sampleid, gene_name), 
                       by = c("method", "sampleid", "gene_name")) %>%
    bind_rows(newa1) %>%
    arrange(method, sampleid, gene_name, lineage)

fig_s4_b <- plot_ranks("HLA-A", allele_newa_df) +
    theme(plot.background = element_rect(fill = "white", color = "white"),
          plot.margin = margin(.5, 1, .5, 1, unit = "cm"))

fig_s4 <- plot_grid(fig_s4_a, fig_s4_b, 
                    ncol = 1, rel_heights = c(1, .7), 
                    labels = c("A)", "B)"))


ggsave("./plots/fig_s4.jpg", fig_s4, width = 6, height = 6)



















## HLA-A isoforms and alleles


ids_t1 <- grep("_t1$", ids, value = TRUE)

salmon_hla_tx <- "./salmon-pers/quant/%s/quant.sf" %>%
    sprintf(ids_t1) %>%
    setNames(sub("_t1$", "", ids_t1)) %>%
    map_df(. %>% 
               read_tsv %>% 
               select(tx_id = Name, counts = NumReads, tpm = TPM) %>%
               filter(grepl("_A\\*", tx_id)) %>%
               separate(tx_id, c("tx_id", "hla"), sep = "_"), 
           .id = "sampleid")

salmon_hla_tx <- salmon_hla_tx %>%
    group_by(sampleid, tx_id) %>%
    mutate(n = ifelse(n() == 1, "1_1", "1")) %>%
    ungroup() %>%
    mutate_at(vars(counts:tpm), ~ifelse(n == "1_1", ./2, .)) %>%
    separate_rows(n, sep = "_") %>%
    select(-n)
    
salmon_hla_tx %>%
    mutate(lineage = sub("^(A\\*\\d+).*$", "\\1", hla)) %>%
    ggplot(aes(reorder_within(lineage, by = tpm, within = tx_id), tpm)) +
    geom_quasirandom(method = "smiley") +
    scale_x_reordered() +
    facet_wrap(~tx_id, ncol = 2, scale = "free") +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = NULL)




















## additional
## 
## 
## 




# New allele z-scores from Veron for HLA-A

nci_allele_a_science <- "/raid/genevol/nci_rnaseq/A_zscores_science.xlsx" %>%
    readxl::read_excel() %>%
    select(lineage = A2d, zall = AZscore)

nci_allele_a_96 <- "/raid/genevol/nci_rnaseq/AZscore_96_66K.xlsx" %>%
    readxl::read_excel() %>%
    select(lineage = Allele, z96 = Azscore) %>%
    mutate(lineage = sub("^(A)(\\d+)$", "\\1*\\2", lineage))

rnaseq_allele_means <- salmon_allele %>%
    filter(gene_name == "HLA-A") %>%
    group_by(lineage) %>%
    summarise(m = mean(rna),
              n = n()) %>%
    ungroup()

new_df <- inner_join(nci_allele_a_science, nci_allele_a_96, by = "lineage") %>%
    inner_join(rnaseq_allele_means, by = "lineage")

cor_qpcr_z <- new_df %>%
    summarise(x = min(zall) * 1.1,
              y = max(z96) * 1.1,
              r = round(cor(z96, zall), 2),
              rho = round(cor(z96, zall, method = "spearman"), 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

z_limits <- range(c(new_df$zall, new_df$z96))

zscores_plot1 <- ggplot(new_df, aes(zall, z96)) +
    geom_abline(linetype = 2) +
    geom_point() +
    geom_point(size = 2) +
    geom_text_repel(aes(label = lineage), size = 2.5) +
    geom_text(data = cor_qpcr_z, aes(x, y, label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3.5) +
    scale_x_continuous(limits = z_limits * 1.1) +
    scale_y_continuous(limits = z_limits * 1.1) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey96"),
          plot.title = element_text(size = 10),
          plot.margin = margin(1, 1, 1, 1, unit = "cm")) +
    labs(x = "z-scores in Science paper", y = "z-score in the 96 donors",
         title = "Correlation of z-scores generated from >500 donors in the Science paper\nvs. z-scores generated in the 96 individuals.")

qpcr_metrics_df <- allele_df %>% 
    filter(gene_name == "HLA-A", method == "qPCR") %>%
    group_by(lineage) %>%
    summarise(lm_mean = mean(rna)) %>%
    ungroup() %>%
    inner_join(new_df, by = "lineage") %>%
    select(lineage, zall, z96, lm_mean) %>%
    pivot_longer(zall:z96, names_to = "z_type", values_to = "z") %>%
    mutate(z_type = recode(z_type, "zall" = "z-scores Science", "z96" = "z-scores 96 donors"))

qpcr_metrics_cor <- qpcr_metrics_df %>%
    group_by(z_type) %>%
    summarise(x = min(lm_mean),
              rho = round(cor(z, lm_mean, method = "spearman"), 2),
              rho_lab = paste("~rho == ", rho)) %>%
    ungroup()

zscores_plot2 <- ggplot(qpcr_metrics_df, aes(lm_mean, z)) +
    geom_point(size = 2) +
    geom_text_repel(aes(label = lineage), size = 2.5) +
    geom_text(data = qpcr_metrics_cor, aes(x = x, y = z_limits[2], label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3.5) +
    facet_wrap(~z_type, nrow = 1) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey96"),
          plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 9),
          plot.margin = margin(1, 1, 1, 1, unit = "cm")) +
    labs(x = "qPCR mean expression (linear model)", y = "qPCR z-score",
         title = "Correlation between qPCR allelic means (linear model) vs. qPCR z-scores.",
         subtitle = "Only including alleles present in at least 5 donors.")






new_long_df <- new_df %>%
    mutate(m_std = GenABEL::rntransform(m)) %>%
    pivot_longer(zall:z96, names_to = "z_type", values_to = "z") %>%
    mutate(z_type = recode(z_type, "zall" = "z-scores Science", "z96" = "z-scores 96 donors"))

cor_new <- new_long_df %>%
    filter(n >= 5) %>%
    group_by(z_type) %>%
    summarise(x = min(m_std),
              rho = round(cor(z, m_std, method = "spearman"), 2),
              rho_lab = paste("~rho == ", rho)) %>%
    ungroup()

zscores_plot3 <- ggplot(filter(new_long_df, n >= 5), aes(m_std, z)) +
    geom_point(size = 2) +
    geom_text_repel(aes(label = lineage), size = 2.5) +
    geom_text(data = cor_new, aes(x = x, y = z_limits[2], label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3.5) +
    facet_wrap(~z_type, nrow = 1) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey96"),
          plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 9),
          plot.margin = margin(1, 1, 1, 1, unit = "cm")) +
    labs(x = "RNA-seq mean expression (std-normal)", y = "qPCR z-score",
         title = "Correlation between RNA-seq allelic means vs. qPCR z-scores.",
         subtitle = "Only including alleles present in at least 5 donors.")

zscores_plot <- 
    plot_grid(plot_grid(zscores_plot1, NULL),
              zscores_plot2, 
              zscores_plot3, 
              ncol = 1) +
    theme(plot.background = element_rect(color = "white", fill = "white"))

ggsave("./plots/zscores.png", zscores_plot, width = 6.5, height = 10)


# new_ranks <- new_df %>%
#     mutate(r_qPCR = rank(-AZscore),
#            r_RNAseq = rank(-m)) %>%
#     pivot_longer(starts_with("r_"), names_to = "method", values_to = "rank") %>%
#     mutate(method = sub("r_", "", method))
# 
# 
# new_p2 <- ggplot(new_ranks, aes(method, rank)) +
#     geom_line(aes(group = lineage), color = "grey") +
#     geom_point() +
#     geom_label(data = filter(new_ranks, method == "qPCR"),
#                aes(label = lineage),
#                nudge_x = -.1) +
#     geom_label(data = filter(new_ranks, method == "RNAseq"),
#                aes(label = lineage),
#                nudge_x = +.1) +
#     scale_y_reverse(breaks = 1:16) +
#     theme_bw() +
#     theme(panel.grid = element_blank()) +
#     labs(x = NULL)
# 
# plot_grid(new_p1, new_p2, nrow = 1)
# ggsave("./plots/newAest.png", width = 8, height = 4)


# Test rankings, compare with Rene et al

rene <- tribble(~lineage, ~r,
                "A*29", 1,
                "A*01", 2, 
                "A*26", 3,
                "A*30", 4,
                "A*03", 5,
                "A*32", 6,
                "A*11", 7,
                "A*25", 8,
                "A*24", 9,
                "A*68", 10,
                "A*23", 11,
                "A*33", 12,
                "A*31", 13,
                "A*02", 14) %>% 
    filter(lineage %in% allele_df$lineage) %>%
    mutate(method = "Rene", r = rank(r)) %>%
    select(method, lineage, r)


allele_df %>%
    filter(gene_name == "HLA-A") %>%
    group_by(method, lineage) %>%
    summarise(m = median(rna)) %>%
    mutate(r = rank(m)) %>%
    ungroup() %>%
    arrange(method, r) %>%
    select(method, lineage, r) %>%
    split(.$method) %>%
    map(~bind_rows(., rene) %>%
            mutate(method = fct_inorder(method)) %>%
            ggplot(aes(method, factor(r))) +
            geom_point() +
            geom_line(aes(group = lineage)) +
            geom_label(data = . %>% filter(method == first(method)),
                       aes(label = lineage), nudge_x = -.15) +
            geom_label(data = . %>% filter(method != first(method)),
                       aes(label = lineage), nudge_x = +.15) +
            theme_minimal() +
            theme(panel.grid = element_blank()) +
            labs(x = NULL, y = NULL)) %>%
    plot_grid(plotlist = ., nrow = 1)









## Fig S5


b2m_temp <- salmon_pers_tpm %>%
    filter(id %in% c("HLA-A", "HLA-B", "HLA-C", "B2M")) %>%
    select(gene_name = id, starts_with("66K")) %>%
    pivot_longer(-gene_name, names_to = "sampleid", values_to = "rna") %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_") %>%
    filter(timepoint == "t1") %>%
    select(sampleid, gene_name, rna)

b2m_df <- b2m_temp %>%
    pivot_wider(names_from = gene_name, values_from = rna) %>%
    pivot_longer(B2M, names_to = "B2M") %>%
    select(sampleid, B2M = value, starts_with("HLA")) %>%
    pivot_longer(starts_with("HLA"), names_to = "gene_name", values_to = "tpm")

b2m_std_df <- b2m_temp %>%
    group_by(gene_name) %>%
    mutate(rna = GenABEL::rntransform(rna)) %>%
    ungroup() %>%
    pivot_wider(names_from = gene_name, values_from = rna) %>%
    pivot_longer(B2M, names_to = "B2M") %>%
    select(sampleid, B2M = value, starts_with("HLA")) %>%
    pivot_longer(starts_with("HLA"), names_to = "gene_name", values_to = "std")

b2m_p1 <- b2m_df %>%
    distinct(sampleid, B2M) %>%
    pivot_longer(B2M, names_to = "gene", values_to = "TPM") %>%
    ggplot(aes(x = gene, y = TPM)) +
    geom_quasirandom(size = 1, method = "smiley") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 9),
          plot.margin = unit(c(.25, 4, .25, 4), "cm")) +
    labs(x = " ", color = "")

# correlation B2m vs HLA
cor_b2m <- b2m_std_df %>%
    group_by(gene_name) %>%
    summarise(x = min(std),
              y = max(B2M),
              r = round(cor(B2M, std), 2),
              rho = round(cor(B2M, std, method = "spearman"), 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

b2m_p2 <- ggplot(b2m_std_df, aes(std, B2M)) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_point(size = .75) +
    geom_text(data = cor_b2m, aes(x, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3) +
    facet_wrap(~gene_name, scales = "free_x") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size = 9)) +
    labs(x = "std norm TPM", y = "std norm TPM (B2M)")

plot_grid(b2m_p1,
          b2m_p2,
          rel_heights = c(1, 1),
          ncol = 1,
          labels = c("A)", "B)"),
          label_size = 11, hjust = 0)

ggsave("./plots/b2m.jpeg", width = 5, height = 3.5)















# Additional

## Fig3 with STD-normalized qPCR values
qpcr_std <- nci_expression %>%
    select(1:3) %>%
    drop_na() %>%
    group_by(gene_name) %>%
    mutate(m_rna = GenABEL::rntransform(m_rna)) %>%
    ungroup()

ggsave("./plots/rnaseq_qpcr_rntrans.jpg", plot1 + labs(title = NULL), height = 2, width = 6)

## Fig 3 showing sample IDs
quant_df2 <- quant_df %>%
    filter(method == "Personalized") %>%
    mutate(sampleid = str_remove(sampleid, "^66K00"))

plot1_sample_ids <- quant_df2 %>%
  ggplot(aes(rnaseq, qPCR)) +
  geom_text(aes(label = sampleid), size = 1, alpha = .5) +
  scale_color_gradient(low = "grey", high = "black") +
  facet_wrap(~gene_name, nrow = 1, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8)) +
  labs(x = "RNA-seq")

ggsave("./plots/rnaseq_qpcr_sampleids.jpg", plot1_sample_ids, height = 2, width = 6)

## qPCR - RNAseq deviations
qpcr_test_plot <- quant_df %>%
  filter(method == "Personalized") %>%
  mutate(sampleid = str_remove(sampleid, "^66K00")) %>%
  mutate(d = qPCR - rnaseq) %>%
  ggplot(aes(x = reorder(sampleid, d, min), y = d)) +
  geom_point(aes(color = gene_name)) +
  scale_color_manual(values = c("HLA-A" = "cornflowerblue", 
			       "HLA-B" = "tomato3", 
			       "HLA-C" = "black")) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
	panel.grid.major.y = element_blank(),
	panel.grid.minor.y = element_blank(),
	axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 1),
	legend.position = "top") +
  labs(x = NULL, y = "qPCR - RNAseq")

ggsave("./plots/rnaseq_qpcr_diffs.jpg", qpcr_test_plot, height = 3, width = 10)

# #Not appropriate to compare loci in qPCR data, but let's see what we get:
qpcr_test_df <- qpcr_std %>%
    pivot_wider(names_from = "gene_name", values_from = "m_rna") %>%
    select(sampleid, `HLA-C`, `HLA-A`, `HLA-B`) %>%
    pivot_longer(-(1:2), names_to = "gene_name") %>%
    drop_na()

cor_qpcr_df <- qpcr_test_df %>%
  group_by(gene_name) %>%
  summarise(y = max(`HLA-C`),
            r = round(cor(`HLA-C`, value), 2),
            rho = round(cor(`HLA-C`, value, method = "spearman"), 2),
            p = cor.test(value, `HLA-C`, method = "pearson", exact = FALSE)$p.value,
            p = format(p, format = "e", digits = 2),
            rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
  ungroup()

out <- ggplot(qpcr_test_df, aes(value, `HLA-C`)) +
    geom_point(size = .75) +
    geom_text(data = cor_qpcr_df, 
            aes(x = -2.5, y + (y *.1), label = rho_lab),
	    hjust = "inward", vjust = "inward", 
            parse = TRUE, size = 3) +
    facet_wrap(~gene_name, switch = "x") +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white"),
	  strip.placement = "outside") +
    labs(x = NULL)

ggsave("./plots/qpcr_test.jpg", out, height = 2, width = 6)

## Fig3 for 11 individuals in batch #2
fig3_11inds <- quant_df %>%
  filter(method == "Personalized") %>%
  mutate(col = sampleid %in% quant_batches$sampleid) %>%
  ggplot(aes(rnaseq, qPCR)) +
  geom_smooth(method = "lm", se = FALSE, color = "grey") +
  geom_point(aes(fill = col), shape = 21, size = 1.5, show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
  facet_wrap(~gene_name, nrow = 1, scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8)) +
  labs(x = "RNA-seq")

ggsave("./plots/rnaseq_qpcr_11inds.jpg", fig3_11inds, height = 2, width = 6)

## Fig 3 including homozygotes only
homoz <- read_tsv("./genos_final.tsv") %>%
    select(1:4) %>%
    group_by(sampleid, gene_name) %>%
    filter(n_distinct(allele) == 1) %>%
    ungroup() %>%
    distinct(sampleid, gene_name)

quant_hom <- quant_df %>%
    filter(method == "Personalized") %>%
    inner_join(homoz) %>%
    select(sampleid, gene_name, qPCR, rnaseq)

cor_hom_df <- quant_hom %>%
    group_by(gene_name) %>%
    summarise(y = max(qPCR),
              r = round(cor(qPCR, rnaseq), 2),
              rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
              p = cor.test(rnaseq, qPCR, method = "pearson", exact = FALSE)$p.value,
              p = format(p, format = "e", digits = 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

ggplot(quant_hom, aes(rnaseq, qPCR)) +
    geom_smooth(method = lm, se = FALSE) +
    geom_point() +
    geom_text(data = cor_hom_df, 
              aes(x = -2.5, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3) +
    facet_wrap(~gene_name, nrow = 1, scales = "free") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8)) +
    labs(x = "RNA-seq")

ggsave("./plots/rnaseq_qpcr_homozygotes.jpg", height = 2, width = 6)

## Fig3 and Fig4 highlighting homozygotes
homoz %>%
    mutate(hom = TRUE) %>%
    left_join(filter(quant_df, method == "Personalized"), .) %>%
    mutate(hom = replace_na(hom, FALSE)) %>%
    select(sampleid, gene_name, hom, qPCR, rnaseq) %>%
    ggplot(aes(rnaseq, qPCR)) +
    geom_point(aes(fill = hom), shape = 21, size = 1.5, stroke = .1) +
    scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
    facet_wrap(~gene_name, nrow = 1, scales = "free") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey96"),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 9)) +
    labs(x = "RNA-seq", fill = "Homozygote?")

ggsave("./plots/rnaseq_qpcr_homozygotes_highlight.jpg", height = 2, width = 6)

homoz %>%
    mutate(hom = TRUE) %>%
    left_join(filter(surface_df, grepl("Personalized", lab)), .) %>%
    mutate(hom = replace_na(hom, FALSE)) %>%
    select(sampleid, gene_name, hom, surface, method, mRNA) %>%
    ggplot(aes(mRNA, surface)) +
    geom_point(aes(fill = hom), shape = 21, stroke = .1, size = 1.5) +
    scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
    facet_wrap(~method, scales = "free_x", strip.position = "bottom") +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.placement = "outside") +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size = 8)) +
    labs(x = NULL, y = "HLA-C surface expression", fill = "Homozygote?")

ggsave("./plots/surface_homozygotes_highlight.jpg", width = 5, height = 2)


## resampling homozygotes
quant_pers_het <- quant_df %>%
    filter(method == "Personalized") %>%
    select(-method, -lab) %>%
    anti_join(homoz)

resamp_df <- quant_pers_het %>% 
    group_by(gene_name) %>%
    nest() %>%
    left_join(count(homoz, gene_name))

resamp_res <- replicate(n = 1000, resamp_df %>%
    mutate(rsp = map2(data, n, ~sample_n(.x, size = .y))) %>%
    select(gene_name, rsp), simplify = FALSE) %>%
    bind_rows(.id = "rep")
    
cor_hom_resamp_df <- resamp_res %>%
    unnest(cols = rsp) %>%
    group_by(rep, gene_name) %>%
    summarise(Pearson = round(cor(qPCR, rnaseq), 2),
              Spearman = round(cor(qPCR, rnaseq, method = "spearman"), 2)) %>%
    ungroup() %>%
    pivot_longer(Pearson:Spearman)

cor_lines_df <- cor_hom_df %>%
    select(gene_name, Pearson = r, Spearman = rho) %>%
    pivot_longer(Pearson:Spearman)

ggplot(cor_hom_resamp_df, aes(x = value, fill = name)) +
    geom_density(size = .25, alpha = .5) +
    geom_vline(data = cor_lines_df, 
               aes(xintercept = value, color = name),
               show.legend = FALSE, linetype = 2) +
    scale_fill_npg() +
    scale_color_npg() +
    facet_wrap(~gene_name) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          plot.background = element_rect(fill = "white"),
          text = element_text(size = 6),
          legend.position = "bottom",
          legend.key.size = unit(0.2, "cm"),
          legend.title = element_text(size = 6),
          plot.caption = element_text(size = 4)) +
    labs(x = "Correlation", fill = "Estimate:",
         title = "Distribution of correlation estimates for heterozygotes",
         subtitle = "Sample sizes: HLA-A (13); HLA-B (8); HLA-C (16)",
         caption = "Vertical lines indicate observed correlation for homozygotes.")
    
ggsave("./plots/correlation_resampling.jpg", width = 4, height = 2.5, dpi = 600)


###

## HLA
hla_batches <- salmon_tx %>%
    filter(sampleid %in% samples_timeps$sampleid) %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    group_by(sampleid, gene_id, gene_name) %>%
    summarise(tpm = sum(tpm)) %>%
    ungroup() %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_") %>%
    mutate(timepoint = toupper(timepoint)) %>%
    pivot_wider(names_from = timepoint, values_from = tpm)

cor_hla_batches <- hla_batches %>%
    group_by(gene_name) %>%
    summarise(min_x = min(T1),
              max_y = max(T2),
              r = round(cor(T2, T1), 2),
              rho = round(cor(T2, T1, method = "spearman"), 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()


### genes
cor_tps_genes <- quant_batches_all_tpm %>%
    group_by(gene_id, gene_name) %>%
    filter(median(t2) > 100) %>%
    summarise(rho = cor(t2, t1, method = "spearman")) %>%
    ungroup()

genes_cor_plot <- cor_tps_genes %>%
    mutate(lab = ifelse(gene_name %in% c("HLA-A", "HLA-B", "HLA-C"), gene_name, "other")) %>%
    ggplot(aes(x = 1, y = rho)) +
    geom_quasirandom(aes(color = lab, size = lab), method = "smiley") +
    scale_color_manual("", values = c("HLA-A" = "blue", "HLA-B" = "tomato", "HLA-C" = "goldenrod3", "other" = "black")) +
    scale_size_manual("", values = c("HLA-A" = 4, "HLA-B" = 4, "HLA-C" = 4, "other" = 1)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          plot.background = element_rect(fill = "white", color = "white")) +
    labs(y = "Spearman's correlation",
         title = "Correlation for all genes with TPM > 100 \nacross 11 individuals in Sample #1 and Sample #2")


cor_grid <- plot_grid(cor_batches_plot + 
                          theme(plot.background = element_rect(fill = "white", color = "white")), 
                      plot_grid(genes_cor_plot, NULL, rel_widths = c(1, 1), nrow = 1) +
                          theme(plot.background = element_rect(fill = "white", color = "white")),
                      ncol = 1,
                      rel_heights = c(1, .5),
                      labels = c("A)", "B)"))

ggsave("./plots/cor_batches.png", cor_grid, width = 8, height = 10)




### HLA
annots <- "../indices/gene_annotation_df.tsv" %>%
    read_tsv() %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    select(gene_name, gene_id)

quant_batches <- read_tsv("./plot_data/peer_residuals.tsv") %>%
    inner_join(annots) %>%
    separate(sampleid, c("sampleid", "batch"), sep = "_") %>%
    mutate(batch = recode(batch, "t1" = "batch 1", "t2" = "batch 2")) %>%
    group_by(sampleid) %>%
    filter(all(c("batch 1", "batch 2") %in% batch)) %>%
    ungroup() %>%
    pivot_wider(names_from = batch, values_from = resid)

quant_batches_std <- quant_batches %>%
    group_by(gene_name) %>%
    mutate_at(vars(4:5), GenABEL::rntransform) %>%
    ungroup()
    
cor_batches <- quant_batches %>%
    group_by(gene_name) %>%
    summarise(min_x = min(`batch 1`),
              max_y = max(`batch 2`),
              r = round(cor(`batch 2`, `batch 1`), 2),
              rho = round(cor(`batch 2`, `batch 1`, method = "spearman"), 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

ggplot(quant_batches, aes(`batch 1`, `batch 2`)) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_point() +
    geom_text(data = cor_batches, 
              aes(x = min_x, y = max_y + .2, label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 2.5) +
    facet_wrap(~gene_name, scales = "free") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size = 8))

ggsave("./plots/timepoints.jpeg", width = 5, height = 2)


ggplot(quant_batches_std, aes(`batch 1`, `batch 2`)) +
    geom_abline() +
    geom_point() +
    facet_wrap(~gene_name) +
    coord_fixed() +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size = 8))



# Remove B2M extreme values
 
rm_b2m_extremes <- b2m_std_df %>%
    filter(abs(B2M) < quantile(abs(B2M), .75)) %>%
    left_join(nci_expression, by = c("sampleid", "gene_name")) %>%
    select(sampleid, gene_name, qPCR = m_rna, RNAseq = std)

cor_b2m_rmExtremes_df <- rm_b2m_extremes %>%
    group_by(gene_name) %>%
    summarise(y = max(qPCR),
              r = round(cor(qPCR, RNAseq), 2),
              rho = round(cor(qPCR, RNAseq, method = "spearman"), 2),
              p = cor.test(RNAseq, qPCR, method = "pearson", exact = FALSE)$p.value,
              p = format(p, format = "e", digits = 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()





# HLA counts normalized by B2M counts
# It does not improve correlation with qPCR
hla_b2m_norm <- read_rds("./plot_data/salmon_pers_b2m.rds") %>%
    left_join(nci_expression, by = c("sampleid", "gene_name")) %>%
    select(sampleid, gene_name, qPCR = m_rna, RNAseq = norm_counts)


cor_hla_b2m_norm_df <- hla_b2m_norm %>%
    group_by(gene_name) %>%
    summarise(y = max(qPCR),
              r = round(cor(qPCR, RNAseq), 2),
              rho = round(cor(qPCR, RNAseq, method = "spearman"), 2),
              p = cor.test(RNAseq, qPCR, method = "pearson", exact = FALSE)$p.value,
              p = format(p, format = "e", digits = 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

hla_b2m_norm %>%
    ggplot(aes(RNAseq, qPCR)) +
    geom_smooth(method = "lm", se = FALSE, color = "skyblue3") +
    geom_point(size = .75) +
    geom_text(data = cor_hla_b2m_norm_df, 
              aes(x = -2.5, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3) +
    facet_wrap(~gene_name, nrow = 1, scales = "free") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 9)) +
    labs(x = "RNA-seq", title = "Personalized reference")


##
quant_b2m_filt <- quant_df_b2m %>%
  filter(between(B2M, -.66, .66))

cor_b2m_filt_df <- quant_b2m_filt %>%
  group_by(gene_name) %>%
  summarise(y = max(qPCR),
            r = round(cor(qPCR, rnaseq), 2),
            rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
            p = cor.test(rnaseq, qPCR, method = "pearson", exact = FALSE)$p.value,
            p = format(p, format = "e", digits = 2),
            rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
  ungroup()

ggplot(quant_b2m_filt, aes(rnaseq, qPCR)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(size = 1.25) +
  scale_color_gradient2() +
  facet_wrap(~gene_name, scales = "free") +
  geom_text(data = cor_b2m_filt_df, aes(x = -2.5, y + (y *.1), label = rho_lab),
            hjust = "inward", vjust = "inward", 
            parse = TRUE,
            size = 3) +
  theme(panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        text = element_text(size = 9),
        plot.margin = unit(c(0, .5, 0, .5), "cm")) +
  labs(x = "RNA-seq", y = "qPCR")

ggsave("./plots/b2m_corr.png", width = 5, height = 2)
ggsave("./plots/b2m_filt.png", width = 5, height = 2)
##
## isoforms

hla_isoforms <- read_tsv("./plot_data/salmon_pers_hla_isoforms.tsv") %>%
  filter(grepl("_t1$", sampleid)) %>%
  filter(gene_name == "HLA-A") %>%
  mutate(sampleid = sub("_t1$", "", sampleid)) %>%
  group_by(sampleid, gene_name, tx_id) %>%
  summarise(tpm = sum(TPM)) %>%
  ungroup() %>%
  group_by(sampleid, gene_name) %>%
  mutate(gene_tpm = sum(tpm)) %>%
  ungroup()

isoforms_df <- left_join(hla_isoforms, 
                         select(nci_expression, sampleid, gene_name, m_rna))

plot_list <- isoforms_df %>%
  split(.$tx_id) %>%
  map(~ggplot(data = ., aes(gene_tpm, m_rna, color = tpm)) +
  geom_point() +
  scale_color_viridis_c() +
  facet_wrap(~tx_id) +
  theme_bw() +
  guides(color = guide_colorbar(barwidth = .5)))

plot_grid(plotlist = plot_list)
ggsave("./plots/hla_a_isoforms.jpeg", width = 10, height = 6)

hla_a_filt <- hla_isoforms %>%
  filter(tx_id != "ENST00000496081.5") %>%
  group_by(sampleid, gene_name) %>%
  summarise(tpm = sum(tpm)) %>%
  ungroup() %>%
  mutate(tpm = scale(tpm)[,1]) %>%
  left_join(select(nci_expression, sampleid, gene_name, m_rna)) %>%
  mutate(sampleid, gene_name, qPCR = m_rna, rnaseq = tpm)

hla_a_filt %>%
  summarise(y = max(qPCR),
            r = round(cor(qPCR, rnaseq), 2),
            rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
            p = cor.test(rnaseq, qPCR, method = "pearson", exact = FALSE)$p.value,
            p = format(p, format = "e", digits = 2),
            rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
  ungroup()



### Transcriptome-wide correlations

# HLA quants from RNA-seq

quant_all <- salmon_pers_tpm %>%
    mutate(gid = sub("\\.\\d+$", "", gid)) %>%
    select(gene_id = gid, gene_name = id, starts_with("66K")) %>%
    pivot_longer(-(gene_id:gene_name), names_to = "sampleid") %>%
    filter(grepl("_t1$", sampleid)) %>%
    mutate(sampleid = sub("_t1$", "", sampleid)) %>%
    group_by(gene_id, gene_name) %>%
    filter(mean(value > 1) > 0.5) %>%
    mutate(value = GenABEL::rntransform(value)) %>%
    ungroup()
    
quant_hla_a <- 
    left_join(filter(quant_all, gene_name == "HLA-A") %>%
                  select(sampleid, hla_gene = gene_name, value),
              quant_all,
              by = "sampleid") 

quant_hla_b <- 
    left_join(filter(quant_all, gene_name == "HLA-B") %>%
                  select(sampleid, hla_gene = gene_name, value),
              quant_all,
              by = "sampleid") 

quant_hla_c <- 
    left_join(filter(quant_all, gene_name == "HLA-C") %>%
                  select(sampleid, hla_gene = gene_name, value),
              quant_all,
              by = "sampleid") 


hla_cors <- 
    bind_rows(quant_hla_a, quant_hla_b, quant_hla_c) %>%    
    group_by(hla_gene, gene_id, gene_name) %>%
    summarise(r = cor(value.x, value.y),
              p = cor.test(value.y, value.x, method = "pearson", exact = FALSE)$p.value) %>%
              
    ungroup() %>%
    filter(hla_gene != gene_name) %>%
    mutate(p_adj = p.adjust(p, method = "BH"))

hla_cors %>%
    filter(p_adj < 0.05) %>%
    group_by(hla_gene) %>%
    top_n(25, abs(r)) %>%
    ungroup() %>%
    ggplot(aes(x = r, y = reorder_within(gene_name, by = r, within = hla_gene))) +
    geom_col() +
    scale_y_reordered() +
    facet_wrap(~hla_gene, nrow = 1, scales = "free") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 7)) +
    labs(x = "Pearson correlation r", y = NULL)

ggsave("./plots/topgenes_hla_correlation.png", width = 6, height = 3)



library(fgsea)

rm_genes <- hla_cors %>%
    distinct(gene_name, gene_id) %>%
    count(gene_name, sort = T) %>%
    filter(n > 1)

gene_lists <- hla_cors %>%
    filter(! gene_name %in% rm_genes$gene_name) %>%
    filter(p_adj < 0.25) %>%
    select(hla_gene, gene_name, r) %>%
    arrange(hla_gene, r) %>%
    split(.$hla_gene) %>%
    map(~select(., -hla_gene)) %>%
    map(deframe)

hallmark <- gmtPathways("./h.all.v7.5.1.symbols.gmt")

gsea_res <- map_df(gene_lists,
                   ~fgsea(pathways = hallmark, stats = .) %>% as_tibble(),
                   .id = "hla_gene")

gsea_res %>%
    mutate(pathway = sub("^HALLMARK_", "", pathway)) %>%
    group_by(hla_gene) %>%
    top_n(25, abs(NES)) %>%
    ungroup() %>%
    ggplot(aes(x = NES, 
               y = reorder_within(pathway, by = NES, within = hla_gene))) +
    geom_col(aes(fill = padj < 0.05), show.legend = FALSE) +
    scale_y_reordered() +
    scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "grey")) +
    facet_wrap(~hla_gene, scales = "free", ncol = 1) + 
    labs(x = "Normalized Enrichment Score",
         y = NULL) + 
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8, margin = margin(r = -2.5)),
          plot.background = element_rect(fill = "white", color = NA),
          plot.caption = element_text(size = 9),
          plot.title = element_text(size = 9)) +
    labs(title = "Top 25 GSEA pathways\n with all correlated genes at FDR = 25%",
         caption = "* Significant enrichment (FDR = 5%) in black.")

ggsave("./plots/gsea.png", width = 4, height = 10)



## protein-coding only?

hla_proteincode_df <- salmon_tx %>%
    filter(timepoint == "Timepoint 1",
           gene_name %in% c("HLA-A", "HLA-B", "HLA-C"),
           tx_type == "protein_coding") %>%
    select(sampleid, tx_id, tx_type, gene_id, gene_name, tpm) %>%
    mutate(sampleid = sub("_t1", "", sampleid)) %>%
    group_by(sampleid, gene_name) %>%
    summarise(tpm = sum(tpm)) %>%
    ungroup() %>%
    left_join(nci_expression) %>%
    select(sampleid, gene_name, rnaseq = tpm, qPCR = m_rna, surface)

hla_proteincode_df %>%
    group_by(gene_name) %>%
    summarise(r = round(cor(qPCR, rnaseq), 2),
              rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
              p = cor.test(rnaseq, qPCR, method = "pearson", exact = FALSE)$p.value,
              p = format(p, format = "e", digits = 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

# coefficient of variation

salmon_pers_hla %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    group_by(gene_name) %>%
    summarise(cv = sd(rna)/mean(rna)) %>%
    ungroup()

nci_expression %>%
    #filter(sampleid %in% unique(salmon_pers_hla$sampleid)) %>%
    drop_na(m_rna) %>%
    group_by(gene_name) %>%
    summarise(cv = sd(m_rna)/mean(m_rna)) %>%
    ungroup()


## size of transcripts?
# expressed_genes <- quant_batches_tpm %>%
#     group_by(gene_id) %>%
#     filter(mean(T1 > 1) > 0.5 | mean(T2 > 1) > 5) %>%
#     ungroup() %>%
#     distinct(gene_id, gene_name)
# 
# len_df <- salmon_tx %>%
#     filter(sampleid %in% samples_timeps$sampleid) %>%
#     inner_join(expressed_genes) %>%
#     group_by(sampleid, timepoint, gene_id, gene_name) %>%
#     summarise(wlen = weighted.mean(len, counts)) %>%
#     ungroup() %>%
#     mutate(timepoint = sub("imepoint ", "", timepoint),
#            sampleid = sub("_t[12]$", "", sampleid)) %>%
#     pivot_wider(names_from = timepoint, values_from = wlen) %>%
#     mutate_at(vars(T1:T2), ~replace_na(., 0)) %>%
#     rename(T1_wlen = T1, T2_wlen = T2)
# 
# 
# len_df %>%
#     filter(T1_wlen > 0 | T2_wlen > 0) %>%
#     ggplot(aes(log10(T1_wlen+1), log10(T2_wlen+1))) +
#     geom_abline(linetype = 2) +
#     geom_point(size = .25) +
#     facet_wrap(~sampleid) +
#     theme_minimal() +
#     theme(panel.grid = element_blank())
# 
# len_df %>%
#     filter(sampleid == "66K00391") %>%
#     filter(T2_wlen == 0 & T1_wlen > 0)


# peer_df <- read_tsv("./plot_data/peer_residuals.tsv") %>%
#     filter(sampleid %in% samples_timeps$sampleid) %>%
#     inner_join(distinct(select_genes, gene_id, gene_name)) %>%
#     separate(sampleid, c("sampleid", "timepoint"), sep = "_") %>%
#     mutate(timepoint = toupper(timepoint)) %>%
#     pivot_wider(names_from = timepoint, values_from = resid)
#     
# cor_tps_peer <- peer_df %>%
#     group_by(sampleid) %>%
#     summarise(rho = round(cor(T2, T1, method = "spearman"), 2),
#               rho_lab = paste("~rho == ", rho)) %>%
#     ungroup()


# TMM
# library(edgeR)
# library(tximport)
# 
# ids_t1 <- grep("t1", ids, value = TRUE)
# 
# files <- "./salmon/quant/%s/quant.sf" %>%
#     sprintf(ids_t1) %>%
#     setNames(ids_t1)
# 
# txdf <- select(tx_annots, TXNAME = tx_id, GENEID = gene_id)
# 
# txi.g <- tximport(files, type="salmon", tx2gene = txdf)
# 
# y <- DGEList(txi.g$counts)
# 
# TMM <- calcNormFactors(y, method = "TMM")
# logCPM <- cpm(TMM, log = TRUE)
# 
# cpm_df <- logCPM %>%
#     as_tibble(rownames = "gene_id") %>%
#     left_join(distinct(tx_annots, gene_id, gene_name)) %>%
#     pivot_longer(starts_with("66K"), names_to = "sampleid", values_to = "logcpm") %>%
#     mutate(sampleid = sub("_t1$", "", sampleid)) %>%
#     select(sampleid, gene_id, gene_name, logcpm)
# 
# tmm_hla <- filter(cpm_df, gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
#     arrange(sampleid, gene_name) %>%
#     left_join(nci_expression) %>%
#     select(sampleid, gene_name, rnaseq = logcpm, qPCR = m_rna)
# 
# cor_tmm <- tmm_hla %>%
#     group_by(gene_name) %>%
#     summarise(x = min(rnaseq),
#               y = max(qPCR), 
#               rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
#               rho_lab = paste("~rho == ", rho)) %>%
#     ungroup()
# 
# ggplot(tmm_hla, aes(rnaseq, qPCR)) +
#     geom_smooth(method = "lm", se = FALSE, color = "skyblue3") +
#     geom_point(size = .75) +
#     geom_text(data = cor_tmm, 
#               aes(x = x, y + (y *.1), label = rho_lab),
#               hjust = "inward", vjust = "inward", 
#               parse = TRUE,
#               size = 2) +
#     facet_wrap(~gene_name, nrow = 1, scales = "free") +
#     theme_bw() +
#     theme(panel.grid.minor = element_blank(),
#           axis.title = element_text(size = 7),
#           axis.text = element_text(size = 6),
#           strip.text = element_text(size = 7),
#           plot.title = element_text(size = 6),
#           plot.margin = margin(.5, .5, 0, .5, unit = "cm")) +
#     labs(x = "RNA-seq", title = "Correlation between qPCR and RNAseq-derived expression levels for Class I HLA genes.")
# 



###

# PCA all samples together, both batches

salmon_all <- salmon_tx %>%
    group_by(sampleid, gene_id, gene_name) %>%
    summarise(tpm = sum(tpm)) %>%
    ungroup()

salmon_all_exp <- salmon_all %>%
    group_by(gene_id, gene_name) %>%
    filter(mean(tpm > 1) > 0.5) %>%
    ungroup()

salmon_matrix <- salmon_all_exp %>%
    select(sampleid, gene_id, tpm) %>%
    pivot_wider(names_from = gene_id, values_from = tpm) %>%
    column_to_rownames("sampleid") %>%
    as.matrix()

pca_all <- prcomp(salmon_matrix, center = TRUE, scale. = TRUE, rank. = 50)

pc_scores_all <- as_tibble(pca_all$x, rownames = "sampleid") %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_") %>%
    mutate(both_tps = sampleid %in% samples_timeps$x) %>%
    select(sampleid, timepoint, both_tps, everything())

ggplot(pc_scores_all, aes(PC1, PC2)) +
    geom_point(aes(fill = timepoint, alpha = both_tps), shape = 21, size = 5) +
    geom_line(data = filter(pc_scores_all, both_tps == TRUE) %>%
                  select(sampleid, timepoint, PC1, PC2),
              aes(group = sampleid), alpha = .5,
              arrow = arrow(length=unit(0.30, "cm"), ends="last", type = "closed")) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = .1)) +
    theme_minimal() +
    theme(panel.grid = element_blank())


pc_loadings_all <- as_tibble(pca_all$rotation, rownames = "gene_id") %>%
    pivot_longer(-gene_id, names_to = "pc") %>%
    mutate(s = ifelse(value > 0, 1, 0)) %>%
    group_by(pc, s) %>%
    slice_max(n = 25, order_by = abs(value)) %>%
    ungroup() %>%
    left_join(distinct(salmon_all_exp, gene_id, gene_name), by = "gene_id") %>%
    select(gene_id, gene_name, pc, value)

pc_loadings_all %>%
    filter(pc %in% c("PC1", "PC2")) %>%
    ggplot(aes(x = value, y = reorder_within(gene_name, by = value, within = pc))) +
    geom_col(aes(fill = value > 0), show.legend = FALSE) +
    scale_y_discrete(labels = function(x) str_remove(x, "_+PC\\d+$")) +
    facet_wrap(~pc, scales = "free") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(size = 10)) +
    labs(x = "PC loadings", y = NULL,
         title = "Top 25 contributions to each direction")

genes_exp <- salmon_all_exp %>%
    distinct(gene_id) %>%
    mutate(gene_id = sub("^([^\\.]+).*$", "\\1", gene_id))
    
# peer



    
peer_matrix <- peer_df %>%
    pivot_wider(names_from = gene_id, values_from = resid) %>%
    column_to_rownames("sampleid") %>%
    as.matrix()

pca_peer <- prcomp(peer_matrix, center = TRUE, scale. = TRUE, rank. = 50)

pc_scores_peer <- as_tibble(pca_peer$x, rownames = "sampleid") %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_") %>%
    mutate(both_tps = sampleid %in% samples_timeps$x) %>%
    select(sampleid, timepoint, both_tps, everything())

ggplot(pc_scores_peer, aes(PC1, PC2)) +
    geom_point(aes(fill = timepoint, alpha = both_tps), shape = 21, size = 5) +
    geom_line(data = filter(pc_scores_peer, both_tps == TRUE) %>%
                  select(sampleid, timepoint, PC1, PC2),
              aes(group = sampleid), alpha = .5,
              arrow = arrow(length=unit(0.30, "cm"), ends="last", type = "closed")) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = .1)) +
    theme_minimal() +
    theme(panel.grid = element_blank())

pc_loadings_peer <- as_tibble(pca_peer$rotation, rownames = "gene_id") %>%
    pivot_longer(-gene_id, names_to = "pc") %>%
    mutate(s = ifelse(value > 0, 1, 0)) %>%
    group_by(pc, s) %>%
    slice_max(n = 10, order_by = abs(value)) %>%
    ungroup() %>%
    left_join(distinct(salmon_all_exp, gene_id, gene_name) %>%
                  mutate(gene_id = sub("^([^\\.]+).*$", "\\1", gene_id)), 
              by = "gene_id") %>%
    select(gene_id, gene_name, pc, value)

pc_loadings_peer %>%
    filter(pc %in% c("PC1", "PC2")) %>%
    ggplot(aes(x = value, y = reorder_within(gene_name, by = value, within = pc))) +
    geom_col(aes(fill = value > 0), show.legend = FALSE) +
    scale_y_discrete(labels = function(x) str_remove(x, "_+PC\\d+$")) +
    facet_wrap(~pc, scales = "free") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(size = 10)) +
    labs(x = "PC loadings", y = NULL,
         title = "Top 10 contributions to each direction")

######
salmon_all_wide <- salmon_all %>%
    filter(sampleid %in% samples_timeps$sampleid) %>%
    mutate(sampleid = sub("_t[12]$", "", sampleid),
           timepoint = sub("imepoint ", "", timepoint)) %>%
    pivot_wider(names_from = timepoint, values_from = tpm) 

tps_means <- salmon_all_wide %>%
    group_by(gene_id, gene_name) %>%
    summarise_at(vars(T1:T2), mean) %>%
    ungroup()

tps_props <- salmon_all_wide %>%
    group_by(sampleid) %>%
    mutate_at(vars(T1:T2), ~./sum(.) * 100) %>%
    ungroup()

ggplot(tps_means, aes(T1, T2)) +
    geom_point()

top_20 <- tps_props %>% 
    mutate(d = (T1 - T2)) %>%
    group_by(gene_id, gene_name) %>%
    summarise(d = mean(d)) %>%
    ungroup() %>%
    arrange(desc(abs(d))) %>%
    slice(1:25)

tps_props %>%
    filter(gene_id %in% top_20$gene_id) %>%
    mutate(gene_id = factor(gene_id, levels = top_20$gene_id)) %>%
    pivot_longer(T1:T2, names_to = "timepoint", values_to = "pct") %>%
    split(.$gene_id) %>%
    map(~ggplot(., aes(sampleid, pct)) +
            geom_col() +
            scale_y_continuous(labels = scales::percent,
                               breaks = scales::pretty_breaks(3)) +
            facet_wrap(~timepoint, ncol = 1) +
            theme_minimal() +
            theme(axis.text.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  plot.title = element_text(size = 8)) +
            labs(x = NULL, y = NULL, title = unique(.$gene_name))) %>%
    plot_grid(plotlist = ., ncol = 5) +
    theme(plot.background = element_rect(color = "white", fill = "white"))


salmon_pers_tx <- "./salmon-pers/quant/%s/quant.sf" %>%
    sprintf(ids) %>%
    setNames(ids) %>%
    map_df(. %>% read_tsv %>% 
               select(tx_id = Name, len = EffectiveLength, counts = NumReads, tpm = TPM), 
           .id = "sampleid") %>%
    mutate(timepoint = ifelse(grepl("t1$", sampleid), "Timepoint 1", "Timepoint 2")) %>%
    mutate(tx_id = ifelse(grepl("_[ABC]\\*", tx_id),
                          sub("^([^_]+).*$", "\\1", tx_id), 
                          tx_id)) %>%
    left_join(tx_annots) %>%
    select(sampleid, timepoint, tx_id, tx_type, len, gene_id, gene_name, counts, tpm)

salmon_pers_hla_tpm <- salmon_pers_tx %>%
    filter(sampleid %in% samples_timeps$sampleid) %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    group_by(sampleid, timepoint, gene_name) %>%
    summarise(tpm = sum(tpm)) %>%
    ungroup() %>%
    mutate(timepoint = sub("imepoint ", "", timepoint),
           sampleid = sub("_t[12]$", "", sampleid)) %>%
    pivot_wider(names_from = timepoint, values_from = tpm)

tps_hla_cor <- salmon_pers_hla_tpm %>%
    group_by(gene_name) %>%
    summarise(x = min(pmin(T1, T2)),
              y = max(pmax(T1, T2)),
              r = round(cor(T2, T1), 2),
              rho = round(cor(T2, T1, method = "spearman"), 2),
              p = cor.test(T1, T2, method = "pearson", exact = FALSE)$p.value,
              p = format(p, format = "e", digits = 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

hla_limits <- salmon_pers_hla_tpm %>%
    pivot_longer(T1:T2, names_to = "dummy") %>%
    group_by(gene_name) %>%
    summarise(r = range(value)) %>%
    ungroup() %>%
    mutate(r2 = r) %>%
    select(gene_name, x = r, y = r2)

tps_plot1 <- ggplot() +
    geom_blank(data = hla_limits, aes(x, y)) +
    facet_wrap(~gene_name, scales = "free") +
    geom_abline(linetype = 2) +
    geom_point(data = salmon_pers_hla_tpm,
               aes(T1, T2), 
               size = 2.5) +
    geom_text_repel(data = salmon_pers_hla_tpm, aes(T1, T2, label = sampleid), 
                    size = 2.5, min.segment.length = 0) +
    geom_text(data = tps_hla_cor, 
              aes(x = x, y, label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 4) +
    scale_color_viridis_c() +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey96"),
          plot.margin = margin(.5, .5, 0, .5, unit = "cm")) +
    labs(x = "Timepoint 1", y = "Timepoint 2",
         title = "Raw TPM estimates")


# peer
peer_df <- read_tsv("./plot_data/peer_residuals.tsv")

hla_ids <- tx_annots %>%
    distinct(gene_id, gene_name) %>% 
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C"))

peer_hla_df <- peer_df %>%
    inner_join(hla_ids) %>%
    inner_join(samples_timeps) %>%
    select(sampleid = x, timepoint, gene_name, resid) %>%
    group_by(gene_name, timepoint) %>%
    mutate(resid = GenABEL::rntransform(resid)) %>%
    ungroup() %>%
    pivot_wider(names_from = timepoint, values_from = resid)

cor_tps_peer_df <- peer_hla_df %>%
    group_by(gene_name) %>%
    summarise(x = min(t1),
              y = max(t2),
              rho = round(cor(t2, t1, method = "spearman"), 2),
              rho_lab = paste("~rho == ", rho)) %>%
    ungroup()

tps_plot2 <- ggplot(peer_hla_df, aes(t1, t2)) +
    geom_abline(linetype = 2) +
    geom_point(size = 2.5) +
    geom_text(data = cor_tps_peer_df, 
              aes(x = x, y * 1.1, label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 4) +
    facet_wrap(~gene_name, nrow = 1, scales = "free") +
    facet_wrap(~gene_name, scales = "free") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey96"),
          plot.margin = margin(.5, .5, 0, .5, unit = "cm")) +
    labs(x = "Timepoint 1", y = "Timepoint 2",
         title = "PEER-corrected estimates for 25 hidden factors, std-normal transformation.")

## DEseq2 normalization

deseq_df <- read_tsv("./plot_data/deseq_counts_separatenorm.tsv") %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    pivot_wider(names_from = timepoint, values_from = scaled_counts)

cor_tps_deseq <- deseq_df %>%
    group_by(gene_name) %>%
    summarise(x = min(t1),
              y = max(t2),
              rho = round(cor(t2, t1, method = "spearman"), 2),
              rho_lab = paste("~rho == ", rho)) %>%
    ungroup()

limits_deseq <- deseq_df %>%
    pivot_longer(t1:t2, names_to = "dummy") %>%
    group_by(gene_name) %>%
    summarise(r = range(value)) %>%
    ungroup() %>%
    mutate(r2 = r) %>%
    select(gene_name, x = r, y = r2)

ggplot() +
    geom_blank(data = limits_deseq, aes(x, y)) +
    facet_wrap(~gene_name, scales = "free") +
    geom_abline(linetype = 2) +
    geom_point(data = deseq_df, aes(t1, t2), size = 2.5) +
    geom_text(data = cor_tps_deseq, 
              aes(x = x, y * 1.1, label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 4) +
    geom_text_repel(data = deseq_df,
                    aes(t1, t2, label = sampleid), 
                    size = 2.5, min.segment.length = 0) +
    facet_wrap(~gene_name, nrow = 1, scales = "free") +
    facet_wrap(~gene_name, scales = "free") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey96"),
          plot.margin = margin(.5, .5, 0, .5, unit = "cm")) +
    labs(x = "Timepoint 1", y = "Timepoint 2",
         title = "DEseq2-scaled counts")



# biotypes

tps_df <- inner_join(salmon_tx, select(samples_timeps, sampleid, x), by = "sampleid") %>%
    select(sampleid = x, timepoint, tx_type, gene_id, gene_name, tpm) %>%
    group_by(sampleid, timepoint, gene_id, gene_name) %>%
    summarise(tpm = sum(tpm)) %>%
    group_by(sampleid, timepoint) %>%
    mutate(p = tpm/sum(tpm)) %>%
    ungroup()

top_10_genes <- tps_df %>%
    group_by(gene_id, gene_name) %>%
    summarise(avg_p = mean(p)) %>%
    ungroup() %>%
    arrange(desc(avg_p)) %>%
    slice(1:10)

tps_final <- tps_df %>%
    mutate(lab = ifelse(gene_id %in% top_10_genes$gene_id, gene_name, "Other"),
           lab = factor(lab, levels = rev(c(top_10_genes$gene_name, "Other"))),
           timepoint = sub("Timepoint ", "", timepoint)) %>%
    group_by(sampleid, timepoint, lab) %>%
    summarise(p = sum(p)) %>%
    ungroup()

cols <- c("grey90", "black", "lightpink", "green4", "red3", "salmon", 
          "lightslateblue", "tomato4", "steelblue1", "steelblue4", "midnightblue")

tps_plot3 <- ggplot(tps_final, aes(timepoint, p, fill = lab)) +
    geom_col(size = 0) +
    scale_fill_manual(values = cols) +
    scale_y_continuous(labels = scales::percent) +
    facet_wrap(~sampleid, nrow = 1) +
    theme_minimal() +
    theme(panel.background = element_rect(color = "white", fill = "white"),
          panel.grid = element_blank(),
          plot.margin = margin(.5, .5, 0, .5, unit = "cm"),
          legend.key.height = unit(.25, "cm"),
          strip.text = element_text(size = 8)) +
    labs(fill = "Gene", y = NULL, 
         title = "Proportion of total transcripts contributed by specific genes.")

tps_biotype <- inner_join(salmon_tx, select(samples_timeps, sampleid, x), by = "sampleid") %>%
    select(sampleid = x, timepoint, tx_type, gene_id, gene_name, tpm) %>%
    group_by(sampleid, timepoint, tx_type) %>%
    summarise(tpm = sum(tpm)) %>%
    group_by(sampleid, timepoint) %>%
    mutate(p = tpm/sum(tpm)) %>%
    ungroup()

top_5_types <- tps_biotype %>%
    group_by(timepoint, tx_type) %>%
    summarise(avg_p = mean(p)) %>%
    ungroup() %>%
    arrange(desc(avg_p)) %>%
    group_by(timepoint) %>%
    top_n(5, avg_p) %>%
    ungroup() %>%
    arrange(timepoint, desc(avg_p)) %>%
    distinct(tx_type) %>%
    pull(tx_type)

tps_biotype_final <- tps_biotype %>%
    mutate(lab = ifelse(tx_type %in% top_5_types, tx_type, "Other"),
           lab = factor(lab, levels = rev(c(top_5_types, "Other"))),
           timepoint = sub("Timepoint ", "", timepoint)) %>%
    group_by(sampleid, timepoint, lab) %>%
    summarise(p = sum(p)) %>%
    ungroup()

tps_plot4 <- ggplot(tps_biotype_final, aes(timepoint, p, fill = lab)) +
    geom_col(size = 0) +
    scale_fill_manual(values = c("grey90", pal_npg()(6)[c(-1, -6)], "black")) +
    scale_y_continuous(labels = scales::percent) +
    facet_wrap(~sampleid, nrow = 1) +
    theme_minimal() +
    theme(panel.background = element_rect(color = "white", fill = "white"),
          plot.margin = margin(.5, .5, 0, .5, unit = "cm"),
          panel.grid = element_blank(),
          strip.text = element_text(size = 8)) +
    labs(y = NULL, fill = "Biotype",
         title = "Proportion of total transcripts in each Gencode biotype.")


tps_out <- plot_grid(tps_plot1, tps_plot2, tps_plot3, tps_plot4, ncol = 1,
                     labels = c("A)", "B)", "C)", "D)")) +
    theme(plot.background = element_rect(color = "white", fill = "white"))

ggsave("./plots/timepoints.png", tps_out, width = 8.5, height = 10)


################################################################################

### edgeR CPMs  
cpm_df <- read_tsv("./plot_data/cpm_edger.tsv") %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C"))


quant_cpm_df <- cpm_df %>%
    select(sampleid, gene_name, cpm) %>%
    left_join(nci_expression, by = c("sampleid", "gene_name")) %>%
    select(sampleid, gene_name, qPCR = m_rna, cpm)

cor_cpm_df <- quant_cpm_df %>%
    group_by(gene_name) %>%
    summarise(x = min(cpm),
              y = max(qPCR),
              r = round(cor(qPCR, cpm), 2),
              rho = round(cor(qPCR, cpm, method = "spearman"), 2),
              p = cor.test(cpm, qPCR, method = "pearson", exact = FALSE)$p.value,
              p = format(p, format = "e", digits = 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

quant_cpm_df %>%
    ggplot(aes(cpm, qPCR)) +
    geom_smooth(method = "lm", se = FALSE, color = "black", size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "white", size = 1) +
    geom_point(size = 1) +
    geom_text(data = cor_cpm_df, 
              aes(x = x, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 2.5) +
    facet_wrap(~gene_name, nrow = 1, scales = "free") +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size = 8, family = "Helvetica"),
          plot.margin = margin(.5, .5, 0, .5, unit = "cm"),
          panel.background = element_rect(fill = "grey96")) +
    labs(x = "RNA-seq", title = "Correlation between qPCR and RNAseq-derived CPMs for Class I HLA genes.")



### HLA fold changes in respect to B2M

b2m_fc <- salmon_pers_tx %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C", "B2M")) %>%
    filter(grepl("_t1$", sampleid)) %>%
    mutate(sampleid = sub("_t1$", "", sampleid)) %>%
    group_by(sampleid, gene_name) %>%
    summarise(counts = sum(counts)) %>%
    ungroup() %>%
    pivot_wider(names_from = gene_name, values_from = counts) %>%
    pivot_longer(starts_with("HLA"), names_to = "gene_name", values_to = "counts") %>%
    mutate(fc = (counts - B2M)/B2M) %>%
    left_join(nci_expression) %>%
    select(sampleid, gene_name, qPCR = m_rna, fc)


cor_fc_df <- b2m_fc %>%
    group_by(gene_name) %>%
    summarise(x = min(fc),
              y = max(qPCR),
              r = round(cor(qPCR, fc), 2),
              rho = round(cor(qPCR, fc, method = "spearman"), 2),
              p = cor.test(fc, qPCR, method = "pearson", exact = FALSE)$p.value,
              p = format(p, format = "e", digits = 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

b2m_fc %>%
    ggplot(aes(fc, qPCR)) +
    geom_smooth(method = "lm", se = FALSE, color = "black", size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "white", size = 1) +
    geom_point(size = 1) +
    geom_text(data = cor_fc_df, 
              aes(x = x, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 2.5) +
    facet_wrap(~gene_name, nrow = 1, scales = "free") +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size = 8, family = "Helvetica"),
          plot.margin = margin(.5, .5, 0, .5, unit = "cm"),
          panel.background = element_rect(fill = "grey96")) +
    labs(x = "RNA-seq", title = "Correlation between qPCR and RNAseq FC in respect to B2M.")


# Shorten isoforms
isos <- 
    "../indices/personalize_transcripts/personalized_transcripts_plusShort.fa" %>%
    Biostrings::readDNAStringSet()

iso_df <- tibble(id = names(isos), len = Biostrings::width(isos)) %>%
    filter(grepl("_A\\*", id)) %>%
    separate(id, c("tx_id", "hla"), sep = "_") %>%
    select(hla, tx_id, len) %>%
    group_by(hla, tx_id) %>%
    mutate(lenclass = case_when(n() == 1 ~ "single",
                                n() == 2 & len == min(len) ~ "short",
                                n() == 2 & len == max(len) ~ "long",
                                TRUE ~ NA_character_)) %>%
    ungroup() %>%
    arrange(hla, tx_id, len)

salmon_pers_short <- "./salmon-pers/quant-shorten/%s/quant.sf" %>%
    sprintf(samples_t1) %>%
    setNames(samples_t1) %>%
    map_df(read_tsv, .id = "sampleid")

salmon_short_hla <- salmon_pers_short %>%
    filter(grepl("_A\\*", Name)) %>%
    separate(Name, c("tx_id", "hla"), sep = "_") %>%
    mutate(sampleid = sub("_t1$", "", sampleid),
           lineage = sub("^(A\\*\\d+).*$", "\\1", hla)) %>%
    select(sampleid, tx_id, hla, lineage, len = Length, 
           elen = EffectiveLength, tpm = TPM, count = NumReads) %>%
    left_join(iso_df)

salmon_short_hla %>%
    filter(lenclass != "single") %>%
    ggplot(aes(x = reorder_within(lineage, by = tpm, within = tx_id, fun = median),
           y = tpm)) +
    geom_quasirandom(aes(color = lenclass), method = "smiley", width = .25) +
    scale_x_reordered() +
    facet_wrap(~tx_id, scale = "free") +
    theme(axis.text.x = element_text(angle = 90))
    



