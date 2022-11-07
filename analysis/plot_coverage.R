library(tidyverse)


cov_plot <- cov_df %>%
    left_join(tmp, by = "pos") %>%
    ggplot(aes(pos, cov, color = is_exon)) +
    geom_line(aes(group = sampleid)) +
    scale_x_continuous(labels = function(x) round(x/1e6, 3)) +
    scale_y_continuous(labels = scales::comma) +
    scale_color_manual(values = c("Yes" = "midnightblue", "No" = "tomato3")) +
    facet_wrap(~gene_name, scales = "free", ncol = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = "Position in chr6 (Mb)", y = "Coverage", color = "Annotated as exon\nin protein-coding\nisoform:")

ggsave("./plots/coverage.png", cov_plot, width = 8, height = 4)








genos %>%
    filter(a01_a11 > 0) %>%
    pull(sampleid) %>%
    write_lines("./plot_data/a01_a11_individuals.txt")
    
    
    
    

mapped_reads <- sprintf("/scratch/vitor/results/%s_n.txt", sampleids) %>%
    setNames(sampleids) %>%
    map_chr(readLines) %>% 
    enframe(name = "sampleid", value = "n") %>%
    mutate(n = as.integer(n))

hla_exons <- read_tsv("./plot_data/hla_exon_coords.tsv") %>%
    mutate(pos = map2(start, end, `:`)) %>%
    select(gene_name, exon, pos) %>%
    unnest(cols = pos)

mean_cov_df <- inner_join(cov_df, hla_exons, by = "pos") %>%
    group_by(sampleid, gene_name) %>%
    summarise(meancov = mean(cov)) %>%
    ungroup() %>%
    inner_join(mapped_reads, by = "sampleid") %>%
    mutate(meancov_per_million = meancov/n * 1e6)

ggplot(mean_cov_df, aes(n, meancov)) +
    geom_point() +
    facet_wrap(~gene_name)


quant_cov_df <- inner_join(mean_cov_df, nci_expression, by = c("sampleid", "gene_name")) %>%
    select(sampleid, gene_name, meancov_per_million, qPCR = m_rna)

cor_cov_df <- quant_cov_df %>%
    group_by(gene_name) %>%
    summarise(y = max(qPCR),
              r = round(cor(qPCR, meancov_per_million), 2),
              rho = round(cor(qPCR, meancov_per_million, method = "spearman"), 2),
              p = cor.test(meancov_per_million, qPCR, method = "pearson", exact = FALSE)$p.value,
              p = format(p, format = "e", digits = 2),
              rho_lab = paste("r ==", r, "*', '~rho == ", rho)) %>%
    ungroup()

quant_cov_df %>%
    ggplot(aes(meancov_per_million, qPCR)) +
    geom_smooth(method = "lm", se = FALSE, color = "black", size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "white", size = 1) +
    geom_point(size = 1) +
    geom_text(data = cor_cov_df, 
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
