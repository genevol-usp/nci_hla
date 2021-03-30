library(tidyverse)

# Annotations

annot <- read_tsv("../../indices/transcript_annotation_df.tsv")

tx_annot <- annot %>%
    select(tx_id, gene_name)

gene_annot <- annot %>%
    distinct(gene_id, gene_name)


# Sample IDs
samples <- readLines("../samples.txt")

# Simulated reads
reads_df <- sprintf("../make_simulation_data/readids/%s.txt", samples) %>%
    setNames(samples) %>%
    map_df(~readLines(.) %>% tibble(readname = .), .id = "sampleid")


# Mapped reads
sam_mapped <- sprintf("./quant/%s_hlaMapped.sam", samples) %>%
    setNames(samples) %>%
    map_df(~read_tsv(., col_names = FALSE) %>%
	   distinct() %>%
	   select(readname = X1, flag = X2, chr = X3, start_m1 = X4, start_m2 = X8,
		  status = X16, targets = X18) %>%
	   mutate_at(vars(status:targets), ~sub("X[SNT]:[ZNi]:", "", .)),
           .id = "sampleid")

count(sam_mapped, flag)

sam_mapped_m1 <- filter(sam_mapped, flag %in% c(83, 99))

sam_mapped_m2 <- filter(sam_mapped, flag %in% c(147, 163)) %>%
    rename(start_m2 = start_m1, start_m1 = start_m2)

sam_mapped_pairs <- 
    full_join(sam_mapped_m1, sam_mapped_m2,
	      by = c("sampleid", "readname", "chr", "start_m1", "start_m2", "status"),
	      suffix = c("_m1", "_m2")) %>%
    mutate(status = ifelse(is.na(targets_m1) | is.na(targets_m2), "Unassigned", status),
	   status = ifelse(!is.na(targets_m1) & !is.na(targets_m2) & targets_m1 != targets_m2, 
			   "Unassigned", status)) %>%
    filter(status == "Assigned") %>%
    select(sampleid, readname, target = targets_m1) %>%
    extract(readname, "tx_id", "[^_]+_(ENST[0-9.]+).*$", remove = FALSE)

mapped_results <- sam_mapped_pairs %>%
    left_join(tx_annot, by = "tx_id") %>%
    rename(true_gene = gene_name) %>%
    left_join(gene_annot, by = c("target" = "gene_id")) %>%
    select(sampleid, readname, true_gene, mapped_gene = gene_name)

mapped_summary <- mapped_results %>% 
    count(sampleid, true_gene, mapped_gene) %>%
    group_by(sampleid, mapped_gene) %>%
    mutate(prop = n/sum(n)) %>%
    group_by(true_gene, mapped_gene) %>%
    summarise(prop = mean(prop)) %>%
    ungroup()

write_rds(mapped_summary, "../plot_data/mapped_reads_summary.rds")


sam_simul <- sprintf("./quant/%s_simulReads.sam", samples) %>%
    setNames(samples) %>%
    map_df(~read_tsv(., col_names = paste0("X", 1:18)) %>%
	   select(readname = X1, flag = X2, chr = X3, start_m1 = X4, start_m2 = X8,
		  status = X16, n_targets = X17, targets = X18) %>%
	   mutate_at(vars(status:targets), ~sub("X[SNT]:[ZNi]:", "", .)) %>%
	   mutate(status = ifelse(grepl("^uT:", status), n_targets, status),
		  n_targets = ifelse(status == "Unassigned_Unmapped", NA, n_targets)),
	   .id = "sampleid")

count(sam_simul, flag)

sam_simul_m1 <- filter(sam_simul, flag %in% c(77, 83, 99, 339, 355))

sam_simul_m2 <- filter(sam_simul, ! flag %in% c(77, 83, 99, 339, 355)) %>%
    mutate(flag1 = case_when(flag == 141 ~ 77,
			     flag == 147 ~ 99,
			     flag == 163 ~ 83,
			     flag == 403 ~ 355,
			     flag == 419 ~ 339)) %>%
    rename(start_m2 = start_m1, start_m1 = start_m2, flag2 = flag)

sam_simul_pairs <- 
    full_join(sam_simul_m1, sam_simul_m2, 
	      by = c("sampleid", "readname", "chr", "start_m1", "start_m2", 
		     "status", "n_targets", "targets", "flag" = "flag1")) %>%
    select(sampleid, readname, status, target = targets) %>%
    extract(readname, "tx_id", "[^_]+_(ENST[0-9.]+).*$", remove = FALSE) %>%
    mutate(rank_status = case_when(status == "Assigned" ~ 0,
				   status == "Unassigned_Secondary" ~ 1,
				   status == "Unassigned_Ambiguity" ~ 2,
				   status == "Unassigned_Unmapped" ~ 3),
	   target = ifelse(rank_status != 0, status, target)) %>%
    group_by(sampleid, readname) %>%
    slice(which.min(rank_status)) %>%
    ungroup() %>%
    select(sampleid, readname, tx_id, target)

simul_results <- sam_simul_pairs %>%
    left_join(tx_annot, by = "tx_id") %>%
    rename(true_gene = gene_name) %>%
    left_join(gene_annot, by = c("target" = "gene_id")) %>%
    mutate(gene_name = ifelse(is.na(gene_name), target, gene_name)) %>%
    select(sampleid, readname, true_gene, mapped_gene = gene_name)

simul_summary <- simul_results %>% 
    count(sampleid, true_gene, mapped_gene) %>%
    group_by(sampleid, true_gene) %>%
    mutate(prop = n/sum(n)) %>%
    group_by(true_gene, mapped_gene) %>%
    summarise(prop = mean(prop)) %>%
    ungroup()


write_rds(simul_summary, "../plot_data/simul_reads_summary.rds")
