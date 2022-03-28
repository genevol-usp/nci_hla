library(tidyverse)
library(scales)
library(ggthemes)
library(ggbeeswarm)
library(cowplot)
library(tidytext)
library(ggrepel)
library(ggsci)
library(corrr)

# Main figs

#### locus overall expression
salmon_pers_tpm <- read_tsv("./salmon-pers/quants.bed.gz") 

salmon_pers_hla <- salmon_pers_tpm %>%
  filter(id %in% c("HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", 
                   "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1")) %>%
  pivot_longer(starts_with("66K"), names_to = "sampleid", values_to = "rna") %>%
  separate(sampleid, c("sampleid", "timepoint"), sep = "_") %>%
  filter(timepoint == "t1") %>%
  select(sampleid, gene_name = id, rna)

salmon_pers_hla %>%
  mutate(gene_name = sub("HLA-", "", gene_name)) %>%
  ggplot(aes(fct_reorder(gene_name, rna, .fun = "median", .desc = TRUE), rna)) +
  geom_quasirandom(size = .5, alpha = .5, method = "smiley") +
  theme_bw() +
  labs(x = NULL, y = "TPM")


salmon_pers_hla %>%
  pivot_wider(names_from = gene_name, values_from = rna) %>%
  select(-sampleid) %>%
  correlate() %>%
  rearrange() %>%
  network_plot(min_cor = .1, curved = FALSE) +
  scale_color_gradient2("r", breaks = seq(0, 1, .2),
                      low = "tomato3", mid = "white", high = "midnightblue",
                      midpoint = .3)

ggsave("./plots/correlations.png", width = 6, height = 3.5)

# 
#   pivot_longer(-term, names_to = "gene_2", values_to = "r") %>%
#   select(gene_1 = term, gene_2, r) %>%
#   mutate(r = ifelse(gene_1 == gene_2, 1, r)) %>%
#   ggplot(aes(gene_1, gene_2)) +
#   geom_tile(aes(fill = r)) +
#   scale_fill_gradient(low = "white", high = "midnightblue") +
#   labs(x = NULL, y = NULL)



ggsave("./plots/salmon_pers_gene.jpg", width = 4, height = 2)

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

salmon_pers_std <- read_tsv("./salmon-pers/quants_std.bed") %>%
  filter(id %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
  select(gene_name = id, starts_with("66K")) %>%
  pivot_longer(-gene_name, names_to = "sampleid", values_to = "rnaseq") %>%
  separate(sampleid, c("sampleid", "timepoint"), sep = "_")

salmon_ref_std <- read_tsv("./salmon/quants_std.bed") %>%
  filter(id %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
  select(gene_name = id, starts_with("66K")) %>%
  pivot_longer(-gene_name, names_to = "sampleid", values_to = "rnaseq")


quant_df <- salmon_pers_std %>%
  filter(timepoint == "t1") %>%
  select(-timepoint) %>%
  bind_rows("Personalized" = ., 
            "Ref transcriptome" = salmon_ref_std, 
            .id = "method") %>%
  left_join(nci_expression, by = c("sampleid", "gene_name")) %>%
  mutate(lab = paste0(gene_name, " (", method, ")")) %>%
  select(sampleid, gene_name, method, lab, qPCR = m_rna, rnaseq, surface)

cor_df <- quant_df %>%
  group_by(lab) %>%
  summarise(y = max(qPCR),
            r = round(cor(qPCR, rnaseq), 2),
            rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
            p = cor.test(rnaseq, qPCR, method = "pearson", exact = FALSE)$p.value,
            p = format(p, format = "e", digits = 2),
            rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
  ungroup()

plot1 <- quant_df %>%
  filter(method == "Personalized") %>%
  ggplot(aes(rnaseq, qPCR)) +
  geom_smooth(method = "lm", se = FALSE, color = "skyblue3") +
  geom_point(size = .75) +
  geom_text(data = cor_df %>% 
              filter(grepl("Person", lab)) %>% 
              mutate(gene_name = sub("^(\\S+).*$", "\\1", lab)), 
            aes(x = -2.5, y + (y *.1), label = rho_lab),
            hjust = "inward", vjust = "inward", 
            parse = TRUE,
            size = 3) +
  facet_wrap(~gene_name, nrow = 1, scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 7)) +
  labs(x = "RNA-seq", title = "Personalized reference")

plot2 <- quant_df %>%
  filter(method != "Personalized") %>%
  ggplot(aes(rnaseq, qPCR)) +
  geom_smooth(method = "lm", se = FALSE, color = "skyblue3") +
  geom_point(size = .75) +
  geom_text(data = cor_df %>% 
              filter(grepl("Ref transcrip", lab)) %>% 
              mutate(gene_name = sub("^(\\S+).*$", "\\1", lab)), 
            aes(x = -2.5, y + (y *.1), label = rho_lab),
            hjust = "inward", vjust = "inward", 
            parse = TRUE,
            size = 3) +
  facet_wrap(~gene_name, nrow = 1, scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 7)) +
  labs(x = "RNA-seq", title = "Reference transcriptome")

plot_grid(plot1, plot2, ncol = 1)


ggsave("./plots/real_data.jpeg", height = 4, width = 6)

# Surface expression
surface_df <- quant_df %>%
  filter(method == "Personalized", 
         gene_name == "HLA-C", 
         !is.na(surface)) %>%
  select(-method) %>%
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
  geom_smooth(method = "lm", se = FALSE, color = "skyblue3") +
  geom_point() +
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
  mutate(lineage = sub("^([^:]+).*$", "\\1", allele)) %>%
  group_by(lineage) %>%
  filter(n() > 5L) %>%
  ungroup() %>%
  select(sampleid, gene_name, lineage, rna = tpm)

nci_allele <- read_rds("./plot_data/nci_allele.rds") %>%
  filter(sampleid %in% unique(salmon_allele$sampleid)) %>%
  select(sampleid, gene_name, lineage = allele, rna) %>%
  group_by(lineage) %>%
  filter(n() > 5L) %>%
  ungroup()

allele_df <- 
  bind_rows("qPCR" = nci_allele, "RNA-seq" = salmon_allele, .id = "method")

allele_ranks <- allele_df %>%
  group_by(method, gene_name, lineage) %>%
  summarise(median = median(rna)) %>%
  group_by(method, gene_name) %>%
  mutate(rk = rank(median)) %>%
  ungroup() %>%
  select(-median) %>%
  pivot_wider(names_from = method, values_from = rk) %>%
  mutate(rk_diff = abs(qPCR - `RNA-seq`)) %>%
  group_by(gene_name) %>%
  mutate(rk_diff = rk_diff/max(rk_diff)) %>%
  ungroup()

# Spearman correlation between expression values
allele_wide_df <- allele_df %>%
  group_by(method, sampleid, gene_name) %>%
  mutate(i = seq_len(n())) %>%
  ungroup() %>%
  pivot_wider(names_from = method, values_from = rna)

cor_df_allelelevel <- allele_wide_df %>%
  group_by(gene_name) %>%
  summarise(x = min(`RNA-seq`),
            y = max(qPCR),
            rho = round(cor(qPCR, `RNA-seq`, method = "spearman"), 2),
            p = cor.test(`RNA-seq`, qPCR, method = "spearman", exact = FALSE)$p.value,
            p = format(p, format = "e", digits = 2),
            rho_lab = paste("~rho ==", rho, "*', '~p == ", p)) %>%
  ungroup()

ggplot(allele_wide_df, aes(`RNA-seq`, qPCR)) +
  geom_smooth(method = "lm", se = FALSE, color = "skyblue3") +
  geom_point(size = .5) +
  geom_text(data = cor_df_allelelevel, 
            aes(x = -2.5, y + (y *.1), label = rho_lab),
            hjust = "inward", vjust = "inward", 
            parse = TRUE,
            size = 3) +
  facet_wrap(~gene_name, scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave("./plots/real_data_allelelevel.jpeg", height = 2, width = 6)


allele_df %>%
  count(method, lineage, rna, sort=T)

# in qPCR inference, all individuals homozygous for
# a given allele will have the same expression levels
allele_df %>%
  filter(lineage == "C*07", method == "qPCR") %>%
  group_by(sampleid) %>%
  filter(n() == 2)

allele_wide_df %>%
  filter(gene_name == "HLA-C") %>%
  arrange(qPCR) %>%
  filter(qPCR < 2) %>%
  print(n=Inf)

#####
#####

allele_pcr <- allele_df %>%
  filter(method == "qPCR") %>%
  ggplot(aes(x = reorder(lineage, rna, FUN = "median"), y = rna)) +
  geom_boxplot(fill = "#6297E770", outlier.color = NA, size = .25) +
  geom_quasirandom(size = .75, alpha = .25, method = "smiley") +
  scale_x_discrete(labels = function(x) sub("^([^*])(\\*)(\\d+)", "\\1\\2\n\\3", x)) +
  facet_wrap(~gene_name, ncol = 1, scales = "free") +
  theme_bw() +
  theme(text = element_text(family = "Times", size = 9),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 9, family = "Times", hjust = 0.5)) +
  labs(x = NULL, y = NULL, title = "qPCR")

allele_rnaseq <- allele_df %>%
  filter(method == "RNA-seq") %>%
  ggplot(aes(x = reorder(lineage, rna, FUN = "median"), y = rna)) +
  geom_boxplot(fill = "#6297E770", outlier.color = NA, size = .25) +
  geom_quasirandom(size = .75, alpha = .25, method = "smiley") +
  scale_x_discrete(labels = function(x) sub("^([^*])(\\*)(\\d+)", "\\1\\2\n\\3", x)) +
  facet_wrap(~gene_name, ncol = 1, scales = "free") +
  theme_bw() +
  theme(text = element_text(family = "Times", size = 9),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 9, family = "Times", hjust = 0.5)) +
  labs(x = NULL, y = NULL, title = "RNA-seq")

alleles_a <- plot_grid(allele_pcr, allele_rnaseq, nrow = 1)
  
alleles_b <- ggplot(allele_ranks, aes(`RNA-seq`, qPCR, color = rk_diff)) +
  geom_abline(alpha = .25, size = 2) +
  geom_point() +
  geom_text_repel(aes(label = lineage), size = 2.5) +
  scale_color_gradient2(low = "blue", mid = "grey25", high = "red",
                        guide = guide_colourbar(direction = "horizontal",
                                                barheight = .25)) +
  facet_wrap(~gene_name, scales="free", nrow = 1) +
  theme_bw() +
  theme(text = element_text(family = "Times", size = 9),
        panel.grid = element_blank(),
        legend.text = element_blank(),
        legend.position = "top",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-5, -5, -5, -5),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(color = expression(Delta~Rank))

plot_grid(alleles_a, NULL, alleles_b, ncol = 1, 
          rel_heights = c(1, .05, .4), 
          labels = c("A)", "B)"),
          label_size = 9, label_fontfamily = "Times")

ggsave("./plots/allele_level_expression.jpg", width = 5, height = 6)






#A*03
dose_df <- quant_df %>%
  filter(method == first(method)) %>%
  select(-method, -surface) %>%
  left_join(alleles_df, by = c("sampleid", "gene_name")) %>%
  pivot_longer(a1:a2, names_to = "i", values_to = "allele") %>%
  select(-i) %>%
  group_by(allele) %>%
  filter(n_distinct(sampleid) >= 15) %>%
  ungroup() %>%
  mutate(dose = 1) %>%
  group_by(sampleid, gene_name, qPCR, rnaseq, allele) %>%
  summarise(dose = sum(dose)) %>%
  ungroup() %>%
  mutate(allele = factor(allele, levels = sort(unique(allele))))


dose_1 <- ggplot(dose_df, aes(rnaseq, qPCR)) +
  geom_line(stat = "smooth", method = "loess", span = 1, se = FALSE, 
            alpha = .5, size = 1, color = "#6297E770") +
  geom_point(aes(color = factor(dose))) +
  scale_color_manual(values = c("0" = "white", "1" = "grey50", "2" = "black")) +
  scale_x_continuous(breaks = pretty_breaks(3)) +
  scale_y_continuous(breaks = pretty_breaks(3)) +
  facet_wrap(~allele, scales="free") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(family = "Times")) +
  labs(x = "RNA-seq", color = "allele\ndose")


alleles_df <- read_rds("./plot_data/salmon_pers_allele.rds") %>%
  filter(grepl("_t1$", sampleid)) %>% 
  mutate(sampleid = sub("_t1$", "", sampleid),
         lineage = sub("^([^:]+).*$", "\\1", allele)) %>%
  select(sampleid, gene_name, lineage) %>%
  group_by(sampleid, gene_name) %>%
  mutate(a = paste0("a", 1:n())) %>%
  ungroup() %>%
  pivot_wider(names_from = a, values_from = lineage)

alleles_quant_df <- quant_df %>%
  filter(method == first(method), gene_name %in% c("HLA-A", "HLA-C")) %>%
  select(-surface, -method) %>%
  left_join(alleles_df, by = c("sampleid", "gene_name")) %>%
  mutate(`A*03` = as.integer(a1 == "A*03") + as.integer(a2 == "A*03"),
         `C*07` = as.integer(a1 == "C*07") + as.integer(a2 == "C*07"))

a03 <- alleles_quant_df %>%
  filter(gene_name == "HLA-A") %>%
  mutate(gene_name = "A*03") %>%
  ggplot(aes(rnaseq, qPCR, fill = as.factor(`A*03`))) +
  geom_point(shape = 21, stroke = .25, size = 2) +
  scale_fill_manual(values = c("0" = "white", "1" = "grey50", "2" = "black")) +
  facet_wrap(~gene_name) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(family = "Times", size = 11)) +
  labs(x = "RNA-seq", fill = "allele\ndose")

c07 <- alleles_quant_df %>%
  filter(gene_name == "HLA-C") %>%
  mutate(gene_name = "C*07") %>%
  ggplot(aes(rnaseq, qPCR, fill = as.factor(`C*07`))) +
  geom_point(shape = 21, stroke = .25, size = 2) +
  scale_fill_manual(values = c("0" = "white", "1" = "grey50", "2" = "black")) +
  facet_wrap(~gene_name) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(family = "Times", size = 11)) +
  labs(x = "RNA-seq", fill = "allele\ndose")


dose_2 <- plot_grid(a03 + theme(legend.position = "none"),
                    c07 + theme(legend.position = "none"),
                    get_legend(a03),
                    nrow = 1, rel_widths = c(1, 1, .25))

plot_grid(dose_1, NULL, dose_2, ncol = 1, 
          rel_heights = c(1, .05, .5), 
          labels = c("A)", "", "B)"),
          label_fontfamily = "Times")
  
ggsave("./plots/alleles_dose.jpg", width = 6, height = 6)

###
alleles_quant_filt <- alleles_quant_df %>%
  filter(!(a1 == "A*03" | a2 == "A*03"),
         !(a1 == "C*07" & a2 == "C*07"))

cor_a3_c7 <- alleles_quant_filt %>%
  mutate(lab = case_when(gene_name == "HLA-A" ~ "A*03",
                         gene_name == "HLA-C" ~ "C*07")) %>%
  group_by(lab) %>%
  summarise(y = max(qPCR),
            r = round(cor(qPCR, rnaseq), 2),
            rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
            p = cor.test(rnaseq, qPCR, method = "pearson", exact = FALSE)$p.value,
            p = format(p, format = "e", digits = 2),
            rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
  ungroup()

a03_filt <- alleles_quant_filt %>%
  filter(gene_name == "HLA-A") %>%
  mutate(gene_name = "A*03") %>%
  ggplot(aes(rnaseq, qPCR)) +
  geom_point(aes(fill = as.factor(`A*03`)), shape = 21, stroke = .25, size = 2) +
  geom_text(data = cor_a3_c7 %>% filter(lab == "A*03"), 
            aes(x = -2.5, y + (y *.1), label = rho_lab),
            hjust = "inward", vjust = "inward", 
            parse = TRUE,
            size = 3) +
  scale_fill_manual(values = c("0" = "white", "1" = "grey50", "2" = "black")) +
  facet_wrap(~gene_name) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(family = "Times", size = 11)) +
  labs(x = "RNA-seq", fill = "allele\ndose")

c07_filt <- alleles_quant_filt %>%
  filter(gene_name == "HLA-C") %>%
  mutate(gene_name = "C*07") %>%
  ggplot(aes(rnaseq, qPCR)) +
  geom_point(aes(fill = as.factor(`C*07`)), shape = 21, stroke = .25, size = 2) +
  geom_text(data = cor_a3_c7 %>% filter(lab == "C*07"), 
            aes(x = -2.5, y + (y *.1), label = rho_lab),
            hjust = "inward", vjust = "inward", 
            parse = TRUE,
            size = 3) +
  scale_fill_manual(values = c("0" = "white", "1" = "grey50", "2" = "black")) +
  facet_wrap(~gene_name) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(family = "Times", size = 11)) +
  labs(x = "RNA-seq", fill = "allele\ndose")

## removing "problematic" alleles only improves correlation slightly


###

# Supplementary figs

# hla-mapper
hlamapper <- read_tsv("./hla-mapper/quants_std.bed") %>% 
  filter(id %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
  select(gene_id = gid, gene_name = id, starts_with("66K")) %>%
  pivot_longer(starts_with("66K"), names_to = "sampleid", values_to = "tpm")

quant_df_hlamapper <- hlamapper %>%
  left_join(nci_expression, by = c("sampleid", "gene_name")) %>%
  select(sampleid, gene_name, rnaseq = tpm, qPCR = m_rna)

cor_df_hlamapper <- quant_df_hlamapper %>%
  group_by(gene_name) %>%
  summarise(y = max(qPCR),
            r = round(cor(qPCR, rnaseq), 2),
            rho = round(cor(qPCR, rnaseq, method = "spearman"), 2),
            rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
  ungroup()

ggplot(quant_df_hlamapper, aes(rnaseq, qPCR)) +
  geom_smooth(method = "lm", se = FALSE, color = "skyblue3") +
  geom_point(size = .75) +
  geom_text(data = cor_df_hlamapper, 
            aes(x = -2.5, y + (y *.1), label = rho_lab),
            hjust = "inward", vjust = "inward", 
            parse = TRUE,
            size = 3) +
  facet_wrap(~gene_name, scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "RNA-seq (hla-mapper)")

ggsave("./plots/real_data_hlamapper.jpeg", height = 2, width = 6)






## Timepoints
quant_batches <- salmon_pers_std %>%
    rename(batch = timepoint) %>%
    mutate(batch = recode(batch, "t1" = "batch 1", "t2" = "batch 2")) %>%
    group_by(sampleid) %>%
    filter(all(c("batch 1", "batch 2") %in% batch)) %>%
    ungroup() %>%
    pivot_wider(names_from = batch, values_from = rnaseq)

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
    theme(panel.grid.minor = element_blank())

ggsave("./plots/timepoints.jpeg", width = 5, height = 2)




# B2M
b2m_df <- salmon_pers_tpm %>%
    filter(id %in% c("HLA-A", "HLA-B", "HLA-C", "B2M")) %>%
    select(gene_name = id, starts_with("66K")) %>%
    pivot_longer(-gene_name, names_to = "sampleid", values_to = "rna") %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_") %>%
    filter(timepoint == "t1") %>%
    select(sampleid, gene_name, rna) %>%
    pivot_wider(names_from = gene_name, values_from = rna) %>%
    pivot_longer(B2M, names_to = "B2M") %>%
    select(sampleid, B2M = value, starts_with("HLA")) %>%
    pivot_longer(starts_with("HLA"), names_to = "gene_name", values_to = "tpm")

b2m_scaled_df <- b2m_df %>%
  distinct(sampleid, B2M) %>%
  mutate(B2M = scale(B2M)[,1])

b2m_p1 <- b2m_scaled_df %>%
  pivot_longer(B2M, names_to = "gene", values_to = "scaled TPM") %>%
  ggplot(aes(x = gene, y = `scaled TPM`, color = `scaled TPM`)) +
  geom_quasirandom(size = 2, method = "smiley") +
  scale_color_gradient2(guide = guide_colourbar(barwidth  = .25)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 9),
        plot.margin = unit(c(.25, .5, .25, 1), "cm"),
        panel.background = element_rect(fill = "grey55")) +
  labs(x = " ", y = "scaled expression", color = "")

quant_df_b2m <- quant_df %>%
  filter(method == "Personalized") %>%
  select(sampleid, gene_name, rnaseq, qPCR) %>%
  left_join(b2m_scaled_df)

cor_personalized_df <- cor_df %>%
  filter(grepl("Person", lab)) %>%
  mutate(gene_name = sub("^(\\S+).*$", "\\1", lab))

b2m_p2 <- ggplot(quant_df_b2m, aes(rnaseq, qPCR)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(aes( color = B2M), size = 1.25) +
  scale_color_gradient2() +
  facet_wrap(~gene_name, scales = "free") +
  geom_text(data = cor_personalized_df, aes(x = -2.5, y + (y *.1), label = rho_lab),
            hjust = "inward", vjust = "inward", 
            parse = TRUE,
            size = 3) +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "grey55"),
        panel.grid = element_blank(),
        legend.position = "none",
        text = element_text(size = 9),
        plot.margin = unit(c(0, .5, 0, .5), "cm")) +
  labs(x = "RNA-seq", y = "qPCR")

# correlation B2m vs HLA
cor_b2m <- b2m_df %>%
    group_by(gene_name) %>%
    summarise(x = min(tpm),
              y = max(B2M),
              r = round(cor(B2M, tpm), 2),
              rho = round(cor(B2M, tpm, method = "spearman"), 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

b2m_p3 <- ggplot(b2m_df, aes(tpm, B2M)) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_point(size = .75) +
    geom_text(data = cor_b2m, aes(x, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3) +
    facet_wrap(~gene_name, scales = "free_x", strip.position = "bottom") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size = 9),
          strip.placement = "outside") +
    labs(x = NULL)

plot_grid(plot_grid(NULL, b2m_p1, NULL, nrow = 1, rel_widths = c(.25, 1, .25)),
          NULL,
          b2m_p2,
          NULL,
          b2m_p3, 
          rel_heights = c(1, .05, 1, .05, 1),
          ncol = 1,
          labels = c("A)", "", "B)", "", "C)"),
          label_size = 11, hjust = 0)
  
ggsave("./plots/b2m.jpeg", width = 5, height = 5)

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


ggplot(cor_pca, aes(PCs, rho, group = gene_name)) +
  geom_line(size = 2, color = "grey25", alpha = .6) +
  geom_point(size = 2, color = "grey25") +
  geom_label_repel(data = filter(cor_pca, PCs == 25),
                   aes(label = gene_name)) +
  scale_y_continuous(breaks = seq(0, .6, .1)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x = "# of PCs for RNA-seq correction",
       y = "Spearman correlation with qPCR estimates")

ggsave("./plots/pca_correction.png", height = 3.5, width = 6)


