library(tidyverse)
library(ggsci)
library(ggthemes)
library(ggridges)
library(cowplot)
library(scales)

# simulated transcripts
ident_df <- read_rds("./plot_data/simulated_transcript_identity.rds") %>%
    mutate(locus = sub("^([^*]+).+$", "HLA-\\1", allele)) %>%
    filter(locus %in% paste0("HLA-", c("A", "B", "C"))) %>%
    select(locus, allele, tx, ident) %>%
    arrange(locus, allele, tx) %>%
    mutate(locus = fct_reorder(locus, ident),
           ident = ident/100)

ident_df %>%
    group_by(locus) %>%
    summarise(median(ident))

ggplot(ident_df, aes(ident, locus)) +
    geom_density_ridges(quantile_lines = TRUE, scale = 1,
                        color = "white", fill = "grey30") +
    theme_minimal() +
    theme(text = element_text(family = "Times"),
          panel.grid.minor.x = element_blank()) +
    scale_x_continuous(labels = scales::percent) +
    labs(x = "Transcript sequence identity with the reference",
         y = NULL)

ggsave("./plots/personalized_identity.jpeg")
    


# expression vs distance
dist_df <- read_rds("./plot_data/weighted_distances.rds")
count_rates <- read_rds("./plot_data/count_rates.rds")

plot_dist_df <- left_join(count_rates, dist_df, by = c("sampleid", "gene_name"))

ggplot(plot_dist_df, aes(wdist, rate_counts, color = method)) +
    geom_hline(yintercept = 1L, linetype = 2, size = 1, color = "grey") +
    geom_line(stat = "smooth", method = "loess", span = 1, se = FALSE, 
              alpha = .4, size = 1, show.legend = FALSE) +
    geom_point(size = 1, alpha = .5) +
    scale_x_continuous(labels = function(x) scales::percent(x, accuracy = 1)) +
    scale_color_manual(values = c("Salmon ref genome" = "black",
                                  "Salmon personalized" = "goldenrod4",
                                  "HLApers" = "tomato4",
                                  "hla-mapper::rna" = "midnightblue")) +
    facet_wrap(~gene_name, scales = "free_x") +
    theme_bw() + 
    theme(text = element_text(family = "Times"), 
          panel.grid = element_blank()) +
    labs(x = "Weighted sequence divergence to the HLA reference allele (%)", 
         y = expression(frac(Estimated~counts, True~counts))) +
    guides(color = guide_legend(override.aes = list(size = 2.5)))
    
ggsave("./plots/accuracy.jpeg", width = 6.5, height = 2)


# Salmon vs hla-mapper
hla_quants <- read_rds("./plot_data/hla_est_counts.rds") %>%
    filter(method %in% c("Salmon personalized", "hla-mapper::rna")) %>%
    pivot_wider(names_from = method, values_from = counts)

cor_hla <- hla_quants %>%
    group_by(gene_name) %>%
    summarise(y = max(`hla-mapper::rna`),
              r = round(cor(`hla-mapper::rna`, `Salmon personalized`), 3),
              rho = round(cor(`hla-mapper::rna`, `Salmon personalized`, 
                              method = "spearman"), 3),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

ggplot(hla_quants, aes(`Salmon personalized`, `hla-mapper::rna`)) +
    geom_point() +
    geom_text(data = cor_hla, aes(x = 0, y + (y *.1), label = rho_lab),
              hjust = "inward", vjust = "inward", 
              parse = TRUE,
              size = 3) +
    facet_wrap(~gene_name, scales = "free") +
    scale_x_continuous(labels = comma, breaks = pretty_breaks(3)) + 
    scale_y_continuous(labels = comma, breaks = pretty_breaks(3)) + 
    theme_bw() +
    theme(text = element_text(size = 8))

ggsave("./plots/salmonpers_vs_hlamapper.jpeg", width = 6, height = 2)

# # Coverage
# 
# coverage_df <- read_rds("./plot_data/coverage.rds") %>%
#     arrange(gene_name, i)
# 
# exons_coords <- coverage_df %>%
#     group_by(gene_name, feature) %>%
#     summarise(start = min(i),
#               end = max(i)) %>%
#     ungroup() %>%
#     arrange(gene_name, start) %>%
#     mutate(manualcolor = ifelse(grepl("intron", feature), "0", "1"))
# 
# ggplot(coverage_df, aes(i, cov, group = sampleid)) +
#     geom_rect(data = exons_coords,
#               aes(xmin = start, xmax = end, ymin = 0, ymax = Inf,
#                   fill = manualcolor), alpha = .2, 
#               inherit.aes = FALSE, show.legend = FALSE) +
#     scale_fill_manual(values = c("0" = NA, "1" = "orange")) +
#     geom_line(alpha = .5) +
#     facet_wrap(~gene_name, scales = "free", ncol = 1) +
#     theme_bw() +
#     theme(panel.grid = element_blank(),
#           axis.ticks.x = element_blank(),
#           axis.text.x = element_blank()) +
#     labs(x = NULL, y = "Read depth")
# 
# ggsave("./plots/coverage.jpeg")

# hla-mapper mappings
simul_reads <- read_rds("./plot_data/simul_reads_summary.rds") %>%
    mutate(mapped_gene = recode(mapped_gene, 
                                "Unassigned_Unmapped" = "Unmap",
                                "Unassigned_Secondary" = "2nd",
                                "Unassigned_Ambiguity" = "Ambig"),
           mapped_gene = factor(mapped_gene, unique(mapped_gene)),
           mapped_gene = fct_relevel(mapped_gene, c("2nd", "Ambig", "Unmap"), after = Inf))
                            

mapped_reads <- read_rds("./plot_data/mapped_reads_summary.rds")

reads_plot_1 <- ggplot(simul_reads, aes(true_gene, mapped_gene)) +
    geom_point(aes(size = prop, color = prop)) +
    scale_color_continuous("", 
                          labels = scales::percent, 
                          breaks = seq(0, 1, by = 0.25),
                          limits = c(0, 1)) +
    scale_size(guide = "none") +
    labs(x = "True Gene", y = "Mapped Gene", 
         title = "Destination of simulated reads")
    
reads_plot_2 <- ggplot(mapped_reads, aes(mapped_gene, true_gene)) +
    geom_point(aes(size = prop, color = prop)) +
    scale_color_continuous("", 
                           labels = scales::percent, 
                           breaks = seq(0, 1, by = 0.25),
                           limits = c(0, 1)) +
    scale_size(guide = "none") +
    labs(x = "Mapped Gene", y = "True Gene", 
         title = "Origin of reads mapped to HLA-A, -B and -C")

plot_grid(reads_plot_1, reads_plot_2, nrow = 2, rel_heights = c(1, .5))

ggsave("./plots/hlamapper_mappings.jpeg", height = 5, width = 6)
