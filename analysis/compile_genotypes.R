library(tidyverse)

samples <- read_tsv("./samples.txt", col_names = FALSE) %>%
    filter(X1 == 1) %>%
    pull(X2)

nci_genos <- read_tsv("/raid/genevol/nci_rnaseq/genotypes_nci_tidy.tsv") %>%
    rename(allele_sanger = allele)

kourami_genos <- sprintf("./kourami/results/%s.result", samples) %>%
    setNames(samples) %>%
    map_df(~read_tsv(., col_names = FALSE), .id = "sampleid") %>%
    mutate(gene_name = sub("^([^*]+).*$", "HLA-\\1", X1)) %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    arrange(sampleid, gene_name, X1) %>%
    group_by(sampleid, gene_name) %>%
    mutate(i = 1:n()) %>%
    ungroup() %>%
    select(sampleid, gene_name, i, allele_kourami = X1)

hlapers_genos <- sprintf("./hlapers/genotypes/%s_t1_genotypes.tsv", samples) %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "sampleid") %>%
    mutate(gene_name = paste0("HLA-", locus),
	   allele = sub("IMGT_", "", allele)) %>%
    arrange(sampleid, gene_name, allele) %>%
    group_by(sampleid, gene_name) %>%
    mutate(i = 1:n()) %>%
    ungroup() %>%
    select(sampleid, gene_name, i, allele_hlapers = allele)

hlamapper <- 
    "/media/storage/genevol/vitor/nci/bam_hlamapper/%s_t1/%s.hla-mapper.log" %>%
    sprintf(samples, samples) %>%
    setNames(samples) %>%
    map_df(. %>% 
	read_lines() %>%
	keep(~grepl("^Possible genotype", .)) %>%
	sub("^Possible genotype for ", "", .) %>%
	tibble(info = .) %>%
	separate(info, c("gene_name", "allele"), sep = ": ") %>%
	filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
	separate_rows(allele, sep = ", "),
	.id = "sampleid") %>%
    arrange(sampleid, gene_name, allele) %>%
    group_by(sampleid, gene_name) %>%
    mutate(i = 1:n()) %>%
    ungroup() %>%
    select(sampleid, gene_name, i, allele_hlamapper = allele)

genos_df <- 
    crossing(sampleid = samples, gene_name = c("HLA-A", "HLA-B", "HLA-C"), i = 1:2) %>%
    left_join(nci_genos) %>%
    left_join(hlapers_genos) %>%
    left_join(kourami_genos) %>%
    left_join(hlamapper)

genos_out <- genos_df %>%
    mutate(hlapers_ok = map2_lgl(allele_sanger, allele_hlapers, ~grepl(.x, .y, fixed = TRUE)),
	   kourami_ok = map2_lgl(allele_sanger, allele_kourami, ~grepl(.x, .y, fixed = TRUE)),
	   hlamapper_ok = map2_lgl(allele_sanger, allele_hlamapper, ~grepl(.x, .y, fixed = TRUE)))

genos_out %>%
    group_by(gene_name) %>%
    summarise_at(vars(ends_with("ok")), mean) %>%
    ungroup()


write_tsv(genos_out, "./genos.tsv")


