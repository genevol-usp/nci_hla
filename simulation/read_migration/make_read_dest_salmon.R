library(tidyverse)
library(furrr)

plan(multisession, workers = 8)

# Considering gene length has the following consequence:
# introns of genes that are also exons of overlapping genes
# will have the same weight for the read origin
gene_annot <- read_tsv("./gene_coords.tsv")
tx_annot <- read_tsv("./transcript_annot.tsv")

flag_table <- 
    tribble(~flag_read, ~flag_mate, ~status, ~pri, ~ori,
	    69, 137, "unmap", FALSE, NA,
	    73, 133, "unmap", FALSE, NA,
	    77, 141, "unmap", FALSE, NA,
	    83, 163, "map", TRUE, "-+",
	    85, 169, "unmap", FALSE, NA,
	    89, 165, "unmap", FALSE, NA,
	    99, 147, "map", TRUE, "+-",
	    325, 393, "unmap", FALSE, NA,
	    329, 389, "unmap", FALSE, NA,
	    339, 419, "map", FALSE, "-+",
	    341, 425, "unmap", FALSE, NA,
	    345, 421, "unmap", FALSE, NA,
	    355, 403, "map", FALSE, "+-")

sample_id <- commandArgs(TRUE)

sam <- "/scratch/vitor/temp/%s_hlareads_salmon.sam" |>
    sprintf(sample_id) |>
    read_tsv(col_names = FALSE) |>
    select(readid = X1, chr = X3, flag = X2, pos = X4, score = X5,
	   cigar = X6, pos_mate = X8, frag_len = X9, seq = X10,
	   nmap = X12, mmindex = X13) |>
    mutate_at(vars(nmap, mmindex), parse_number)

first_reads <- sam |> 
    select(readid, flag, chr, pos, score, cigar, pos_mate, mmindex) |>
    inner_join(flag_table, join_by(flag == flag_read)) |>
    filter(status == "map") |>
    select(-status)

second_reads <- sam |> 
    select(readid, flag, chr, pos, score, cigar, pos_mate, mmindex) |>
    inner_join(flag_table, by = join_by(flag == flag_mate)) |>
    filter(status == "map") |>
    rename("flag_mate" = "flag") |>
    select(-status, -flag_read)

alig_df <- 
    left_join(first_reads, second_reads,
	      by = c("readid", "chr", "score", "pos_mate" = "pos", "pos" = "pos_mate",
		     "flag_mate", "mmindex", "pri", "ori"),
	      suffix = c("_read", "_mate")) |>
    rowid_to_column("alig")
    
destinations_df <- alig_df |>
    select(alig, readid, chr) |>
    left_join(tx_annot, by = c("chr" = "transcript_id")) |>
    group_by(readid) |>
    mutate(pct = 1/n()) |>
    group_by(readid, gene_name) |>
    summarise(pct = sum(pct)) |>
    ungroup() |>
    mutate(pos = NA)

# unmap
flags_unmap <- flag_table |>
    filter(status == "unmap") |>
    pivot_longer(starts_with("flag")) |>
    pull(value)

unmap_df <- sam |>
    filter(flag %in% flags_unmap) |>
    distinct(readid) |>
    mutate(gene_name = "unmap", pct = 1, pos = 1)

out <- bind_rows(destinations_df, unmap_df)

write_tsv(out, sprintf("./read_destinations/%s_salmon.tsv", sample_id))
