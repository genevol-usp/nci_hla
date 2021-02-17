library(tidyverse)
library(ggthemes)

# simulated transcripts
ident_df <- read_rds("./plot_data/simulated_transcript_identity.rds") %>%
    mutate(locus = sub("^([^*]+).+$", "HLA-\\1", allele)) %>%
    select(locus, allele, tx, ident) %>%
    arrange(locus, allele, tx) %>%
    mutate(locus = fct_reorder(locus, ident),
           ident = ident/100)

ident_df %>%
    group_by(locus) %>%
    summarise(median(ident))

ggplot(ident_df, aes(ident, locus)) +
    ggridges::geom_density_ridges(quantile_lines = TRUE,
                                  color = "white",
                                  fill = "grey30",
                                  scale = 1) +
    theme_minimal() +
    theme(text = element_text(family = "Times"),
          panel.grid.minor.x = element_blank()) +
    scale_x_continuous(labels = scales::percent) +
    coord_cartesian(xlim = c(.9, 1)) +
    labs(x = "Transcript sequence identity with the reference",
         y = NULL)

ggsave("./plots/personalized_identity.jpeg", width = 5)
    


# MARDs
mard_genome <- read_rds("./plot_data/mard_genomewide.rds")
mard_hla <- read_rds("./plot_data/mard_hla.rds")

mard_df <- bind_rows(genome = mard_genome, hla = mard_hla, .id = "subset")

ggplot(mard_df, aes(mard, subset, fill = method, color = method)) +
    ggridges::geom_density_ridges(scale = 1, alpha = .5) +
    theme_bw() +
    scale_fill_colorblind() +
    scale_color_colorblind() +
    labs(x = "MARD", y = NULL)
