library(tidyverse)
library(scales)
library(ggthemes)
library(ggbeeswarm)
library(cowplot)


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


### Salmon pers
salmon_pers <- read_tsv("./salmon-pers/quants_std.bed") %>%
    filter(id %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    select(gene_name = id, starts_with("66K")) %>%
    pivot_longer(-gene_name, names_to = "sampleid", values_to = "rnaseq") %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_")

nci_expression <- 
    "/raid/genevol/nci_rnaseq/phase1/HLA-A, -B, -C expression levels for Pat.xlsx" %>%
    readxl::read_excel() %>%
    select(sampleid = 3, 5:8) %>%
    mutate(across(-sampleid, as.numeric)) %>%
    pivot_longer(-sampleid, names_to = "tmp") %>%
    separate(tmp, c("gene_name", "experiment"), sep = " ") %>%
    pivot_wider(names_from = experiment, values_from = value)

quant_df <- salmon_pers %>%
    filter(timepoint == "t1") %>%
    select(-timepoint) %>%
    left_join(nci_expression, by = c("sampleid", "gene_name")) %>%
    select(sampleid, gene_name, qPCR = mRNA, rnaseq, surface = Surface)

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

