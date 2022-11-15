library(tidyverse)

samples <- read_lines("../sample_ids_t1.txt")

juncs <- sprintf("/scratch/vitor/splicing/%s.junc", samples)

if (all(file.exists(juncs))) {
    write_lines(juncs, "./juncfiles.txt")
} else {
    stop("missing junction files")
}

genos <- read_tsv("../genos_final.tsv") %>%
    select(sampleid, gene_name, i, allele) %>%
    mutate(lineage = sub("^([ABC]\\*\\d+).*$", "\\1", allele))

grp_data <- genos %>%
    filter(gene_name == "HLA-A") %>%
    group_by(sampleid) %>%
    summarise(grp = ifelse(any(lineage %in% c("A*01", "A*11")), "a1", "ctrl")) %>%
    ungroup() %>%
    mutate(grp = factor(grp, levels = c("ctrl", "a1"))) %>%
    arrange(grp, sampleid)

write_tsv(grp_data, "./group_file.txt", col_names = FALSE)

