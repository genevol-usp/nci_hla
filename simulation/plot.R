library(tidyverse)
library(ggsci)
library(ggthemes)
library(ggridges)
library(cowplot)
library(scales)
library(tidytext)
library(ggforce)

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

plot_dist_df <- left_join(count_rates, dist_df, by = c("sampleid", "gene_name")) %>%
    filter(method != "Salmon hla-mapper")

fig1 <- ggplot(plot_dist_df, aes(wdist, rate_counts)) +
    geom_line(aes(group = method), 
              stat = "smooth", method = "loess", span = 1, se = FALSE, 
              size = 1.5, show.legend = FALSE) +
    geom_line(aes(color = method), 
              stat = "smooth", method = "loess", span = 1, se = FALSE, 
              size = 1, show.legend = FALSE) +
    geom_point(aes(fill = method), size = 2, shape = 21) +
    scale_x_continuous(breaks = pretty_breaks(3),
                       labels = function(x) scales::percent(x, accuracy = 1)) +
    scale_color_manual(labels = c("Salmon ref genome" = "Ref transcriptome",
                                 "Salmon personalized" = "Personalized"),
                      values = c("white", "black")) +
    scale_fill_manual(labels = c("Salmon ref genome" = "Ref transcriptome",
                                  "Salmon personalized" = "Personalized"),
                       values = c("white", "black")) +
    facet_wrap(~gene_name, scales = "free_x") +
    theme(text = element_text(size = 8, family = "Arial"), 
          legend.text = element_text(size = 8, family = "Arial"),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey96"),
          legend.position = "top") +
    labs(x = "Weighted sequence divergence to the HLA allele in the reference genome (%)", 
         y = expression(frac(Estimated~counts, True~counts)),
         fill = "Method:") +
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1)))
    
ggsave("./plots/fig1.jpg", fig1, width = 5, height = 2.5)

# Origin-destination
od_star_data <- read_tsv("./read_migration/origin_dest_results.tsv")

od_star_avg <- od_star_data |>
    complete(sampleid, origin, dest, fill = list(w_by_simul = 0, w_by_mapped = 0)) |>
    group_by(origin, dest) |>
    summarise_at(vars(starts_with("w_by")), mean) |>
    ungroup()

dat_star_simul_gg <- od_star_avg |>
    filter(origin %in% c("HLA-A", "HLA-B", "HLA-C")) |>
    select(origin, dest, w = w_by_simul) |>
    mutate(dest = ifelse(w >= 0.01 | dest == "unmap", dest, "other")) |>
    group_by(origin, dest) |>
    summarise(w = sum(w)) |>
    ungroup() |>
    select(origin, dest, Freq = w) |>
    gather_set_data(1:2) %>%
    arrange(x, origin, desc(Freq))

dat_star_mapped_gg <- od_star_avg |>
    filter(dest %in% c("HLA-A", "HLA-B", "HLA-C")) |>
    select(origin, dest, w = w_by_mapped) |>
    mutate(origin = ifelse(w >= 0.01, origin, "other")) |>
    group_by(origin, dest) |>
    summarise(w = sum(w)) |>
    ungroup() |>
    select(origin, dest, Freq = w) |>
    gather_set_data(1:2) %>%
    arrange(x, origin, desc(Freq))

od_salmon <- read_tsv("./read_migration/origin_dest_results_salmon.tsv")

od_salmon_avg <- od_salmon |>
    complete(sampleid, origin, dest, fill = list(w_by_simul = 0, w_by_mapped = 0)) |>
    group_by(origin, dest) |>
    summarise_at(vars(starts_with("w_by")), mean) |>
    ungroup()

dat_salmon_simul_gg <- od_salmon_avg |>
    filter(origin %in% c("HLA-A", "HLA-B", "HLA-C")) |>
    select(origin, dest, w = w_by_simul) |>
    mutate(dest = ifelse(w >= 0.01 | dest == "unmap", dest, "other")) |>
    group_by(origin, dest) |>
    summarise(w = sum(w)) |>
    ungroup() |>
    select(origin, dest, Freq = w) |>
    gather_set_data(1:2) %>%
    arrange(x, origin, desc(Freq))

dat_salmon_mapped_gg <- od_salmon_avg |>
    filter(dest %in% c("HLA-A", "HLA-B", "HLA-C")) |>
    select(origin, dest, w = w_by_mapped) |>
    mutate(origin = ifelse(w >= 0.01, origin, "other")) |>
    group_by(origin, dest) |>
    summarise(w = sum(w)) |>
    ungroup() |>
    select(origin, dest, Freq = w) |>
    gather_set_data(1:2) %>%
    arrange(x, origin, desc(Freq))

mycols <- 
    c(dat_star_simul_gg$y, dat_star_mapped_gg$y, dat_salmon_simul_gg$y, dat_salmon_mapped_gg$y) |>
    unique() |>
    tibble(dest = _) |>
    mutate(color = case_when(dest == "HLA-A" ~ "black",
			     dest == "HLA-B" ~ "cornflowerblue",
			     dest == "HLA-C" ~ "goldenrod3",
			     TRUE ~ "grey")) |>
    arrange(dest) |>
    deframe()


plot_ggforce <- function(dat) {
    ggplot(dat, aes(x = x, id = id, split = y, value = Freq)) +
	geom_parallel_sets(aes(fill = dest), alpha = .75, axis.width = 0.2,
			   n = 100, strength = 0.5) +
	geom_parallel_sets_axes(axis.width = 0.25, fill = "white",
				color = NA, size = 0.15) +
	geom_parallel_sets_labels(color = "gray35", size = 4, angle = 0, fontface = "bold") +
	scale_fill_manual(values = mycols) +
	scale_color_manual(values = mycols) +
	theme( 
	      panel.background = element_rect(fill = "white", color = "white"),
	      plot.background = element_rect(color = "white", fill = "white"),
	      legend.position = "none",
	      panel.grid.major = element_blank(),
	      panel.grid.minor = element_blank(),
	      strip.background = element_rect(color = "white", fill = NA, size = NA),
	      axis.text = element_blank(),
	      axis.title = element_blank(),
	      axis.ticks = element_blank(),
	      plot.title = element_text(hjust = 0.5))
}

od_out <- 
    plot_grid(
	      plot_ggforce(dat_star_simul_gg) + labs(title = "STAR (weighted by total reads simulated from origin)"),
	      plot_ggforce(dat_star_mapped_gg) + labs(title = "STAR (weighted by total reads mapped to destination)"),
	      plot_ggforce(dat_salmon_simul_gg) + labs(title = "Salmon (weighted by total reads simulated from origin)"),
	      plot_ggforce(dat_salmon_mapped_gg) + labs(title = "Salmon (weighted by total reads mapped to destination)"),
	      ncol = 2)	  

ggsave("./plots/readmigration.jpg", od_out, dpi = 600, height = 10, width = 10)






# Divergent individuals at HLA-B

diverg_b <- dist_df |>
    filter(gene_name == "HLA-B") |>
    filter(wdist > median(wdist)) |>
    arrange(wdist)

od_salmon_div <- filter(od_salmon, sampleid %in% diverg_b$sampleid)

od_salmon_div_avg <- od_salmon_div |>
    complete(sampleid, origin, dest, fill = list(w_by_simul = 0, w_by_mapped = 0)) |>
    group_by(origin, dest) |>
    summarise_at(vars(starts_with("w_by")), mean) |>
    ungroup()

dat_salmon_div_simul_gg <- od_salmon_div_avg |>
    filter(origin %in% c("HLA-A", "HLA-B", "HLA-C")) |>
    select(origin, dest, w = w_by_simul) |>
    mutate(dest = ifelse(w >= 0.01 | dest == "unmap", dest, "other")) |>
    group_by(origin, dest) |>
    summarise(w = sum(w)) |>
    ungroup() |>
    select(origin, dest, Freq = w) |>
    gather_set_data(1:2) %>%
    arrange(x, origin, desc(Freq))

dat_salmon_div_mapped_gg <- od_salmon_div_avg |>
    filter(dest %in% c("HLA-A", "HLA-B", "HLA-C")) |>
    select(origin, dest, w = w_by_mapped) |>
    mutate(origin = ifelse(w >= 0.01, origin, "other")) |>
    group_by(origin, dest) |>
    summarise(w = sum(w)) |>
    ungroup() |>
    select(origin, dest, Freq = w) |>
    gather_set_data(1:2) %>%
    arrange(x, origin, desc(Freq))

    plot_grid(
	      plot_ggforce(dat_salmon_div_simul_gg) + labs(title = "Salmon (weighted by total reads simulated from origin)"),
	      plot_ggforce(dat_salmon_div_mapped_gg) + labs(title = "Salmon (weighted by total reads mapped to destination)"),
	      ncol = 2)	  




# Salmon vs hla-mapper
# hla_quants <- read_rds("./plot_data/hla_est_counts.rds") %>%
#     filter(method %in% c("Salmon personalized", "hla-mapper::rna")) %>%
#     pivot_wider(names_from = method, values_from = counts)
# 
# cor_hla <- hla_quants %>%
#     group_by(gene_name) %>%
#     summarise(y = max(`hla-mapper::rna`),
#               r = round(cor(`hla-mapper::rna`, `Salmon personalized`), 3),
#               rho = round(cor(`hla-mapper::rna`, `Salmon personalized`, 
#                               method = "spearman"), 3),
#               rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
#     ungroup()
# 
# ggplot(hla_quants, aes(`Salmon personalized`, `hla-mapper::rna`)) +
#     geom_point() +
#     geom_text(data = cor_hla, aes(x = 0, y + (y *.1), label = rho_lab),
#               hjust = "inward", vjust = "inward", 
#               parse = TRUE,
#               size = 3) +
#     facet_wrap(~gene_name, scales = "free") +
#     scale_x_continuous(labels = comma, breaks = pretty_breaks(3)) + 
#     scale_y_continuous(labels = comma, breaks = pretty_breaks(3)) + 
#     theme_bw() +
#     theme(text = element_text(size = 8))
# 
# ggsave("./plots/salmonpers_vs_hlamapper.jpeg", width = 6, height = 2)

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

# Other pipelines
plot_dist_alt <- 
    left_join(count_rates, dist_df, by = c("sampleid", "gene_name")) %>%
    filter(!grepl("Salmon", method))

alt_pipes_1 <- ggplot(plot_dist_alt, aes(wdist, rate_counts, color = method)) +
    geom_hline(yintercept = 1L, linetype = 2, size = 1, color = "grey") +
    geom_line(stat = "smooth", method = "loess", span = 1, se = FALSE, 
              alpha = .4, size = 1, show.legend = FALSE) +
    geom_point(size = 1, alpha = .5) +
    scale_x_continuous(labels = function(x) scales::percent(x, accuracy = 1)) +
    scale_color_manual(values = c("HLApers" = "grey25",
                                  "hla-mapper::rna" = "royalblue",
                                  "hla-mapper:rna_corrected" = "goldenrod3")) +
    facet_wrap(~gene_name, scales = "free_x") +
    theme_bw() + 
    theme(text = element_text(family = "Times"), 
          panel.grid = element_blank(),
          legend.position = "top") +
    labs(x = "Weighted sequence divergence to the HLA reference allele (%)", 
         y = expression(frac(Estimated~counts, True~counts))) +
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1)))




simul_reads <- read_rds("./plot_data/simul_reads_summary.rds") %>%
    mutate(mapped_gene = recode(mapped_gene, 
                                "Unassigned_Unmapped" = "Unmap",
                                "Unassigned_Secondary" = "2nd",
                                "Unassigned_Ambiguity" = "Ambig"),
           mapped_gene = factor(mapped_gene, unique(mapped_gene)),
           mapped_gene = fct_relevel(mapped_gene, c("2nd", "Ambig", "Unmap"), after = Inf))
                            

mapped_reads <- read_rds("./plot_data/mapped_reads_summary.rds")

alt_pipes_2 <- ggplot(mapped_reads, aes(mapped_gene, true_gene)) +
    geom_point(aes(size = prop)) +
    scale_color_continuous("", 
                           labels = scales::percent, 
                           breaks = seq(0, 1, by = 0.25),
                           limits = c(0, 1)) +
    scale_size(guide = "none") +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm")) +
    labs(x = "Mapped Gene", y = "True Gene")


alt_pipes_3 <- simul_reads %>%
    filter(prop > 0.0001) %>%
    ggplot(aes(x = reorder_within(mapped_gene, desc(prop), true_gene), 
               y = prop)) +
    geom_col() +
    facet_wrap(~true_gene, ncol = 1, scales = "free_x") +
    scale_x_discrete(labels = function(x) str_extract(x, "^([^_]+)")) + 
    scale_y_continuous(labels = percent, breaks = c(0, .5, 1)) +
    theme(axis.text.x = element_text(size = 7)) +
    labs(x = NULL, y = " ")


subread_summary <- read_rds("./plot_data/summary_subread_mappings.rds") %>%
    mutate(mapped_gene = recode(mapped_gene, 
                                "Unassigned_Unmapped" = "Unmap",
                                "Unassigned_Secondary" = "2nd",
                                "Unassigned_Ambiguity" = "Ambig",
                                "Not present in SAM" = "Not in SAM"))

mapped_levels <- subread_summary %>% 
    group_by(mapped_gene) %>% 
    summarise(prop = mean(prop)) %>%
    ungroup() %>% 
    arrange(desc(prop)) %>%
    pull(mapped_gene)

alt_pipes_4 <- subread_summary %>%
    mutate(mapped_gene = factor(mapped_gene, levels = mapped_levels)) %>%
    ggplot(aes(i, prop, fill = mapped_gene)) +
    geom_col(width = 1, position = "fill") +
    facet_wrap(~gene_name, ncol = 1, scales = "free_x") +
    scale_y_continuous(labels = percent) +
    scale_fill_npg() +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.text = element_blank(), 
          strip.background = element_rect(fill = "grey70")) +
    labs(x = NULL, y = NULL, fill = NULL)


plot_grid(alt_pipes_1, 
          alt_pipes_2, 
          plot_grid(alt_pipes_3, alt_pipes_4, nrow = 1, rel_widths = c(.9, 1)),
          ncol = 1, 
          rel_heights = c(.5, .5, 1),
          labels = c("A)", "B)", "C)"))

ggsave("./plots/alternative_pipelines.jpeg", height = 8, width = 6)







