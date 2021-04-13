library(tidyverse)
library(scales)
library(ggthemes)

salmon <- read_rds("./plot_data/salmon_quants.rds")

tpm_proportions <- salmon %>%
    arrange(sampleid, -tpm) %>%
    group_by(sampleid) %>%
    mutate(prop = tpm/sum(tpm)) %>%
    ungroup() %>%
    mutate(labl = ifelse(prop >= 0.05, gene_name, "Other")) %>%
    group_by(sampleid, gene = labl) %>%
    summarise(prop = sum(prop)) %>%
    ungroup()

gene_levels <- tpm_proportions %>%
    group_by(gene) %>%
    summarise(prop = mean(prop)) %>%
    ungroup() %>%
    arrange(-prop) %>%
    pull(gene)

tpm_proportions <- tpm_proportions %>%
    mutate(gene = factor(gene, levels = gene_levels))

mycols <- c("white", "#046C9A", "#D69C4E", "#ABDDDE",
            "tomato3", "grey60", "grey35", "#ECCBAE")

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
