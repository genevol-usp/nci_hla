library(tidyverse)

get_positions <- function(position, cigar) {
    
    res <- str_extract_all(cigar, "(\\d+[A-Z])") |>
	unlist()

    reps <- as.integer(str_extract(res, "(\\d+)"))
    lets <- str_extract(res, "[A-Z]")

    pattern <- map2(lets, reps, ~rep(.x, .y)) |>
	unlist()

    positions <- position:(position + length(pattern) - 1)

    tibble(position = positions[pattern == "M"])
}

# Considering gene length has the following consequence:
# introns of genes that are also exons of overlapping genes
# will have the same weight for the read origin
gene_annot <- read_tsv("./gene_coords.tsv")

flag_table <- 
    tribble(~flag_read, ~flag_mate, ~status, ~pri, ~ori,
	    77, 141, "unmap", FALSE, NA,
	    83, 163, "map", TRUE, "-+",
	    99, 147, "map", TRUE, "+-",
	    339, 419, "map", FALSE, "-+",
	    355, 403, "map", FALSE, "+-")

sample_id <- commandArgs(TRUE)

sam <- "/scratch/vitor/temp/%s_hlareads.sam" |>
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

first_res <- alig_df |> 
    select(alig, readid, chr, pos, cigar = cigar_read) |>
    mutate(positions = map2(pos, cigar, get_positions)) |>
    select(alig, readid, chr, positions)

second_res <- alig_df |> 
    select(alig, readid, chr, pos = pos_mate, cigar = cigar_mate) |>
    mutate(positions = map2(pos, cigar, get_positions)) |>
    select(alig, readid, chr, positions)


# 7 minutes and 20GB of RAM
# because grouping step is too slow
destinations_read <- first_res |>
    mutate(readlen = map_int(positions, nrow)) |>
    unnest(cols = c(positions)) |>
    left_join(gene_annot, join_by(chr, between(position, start, end))) |>
    add_count(alig, readid, chr, position) |>
    mutate(w = (1/n)/readlen) |>
    group_by(alig, readid, gene_name) |>
    summarise(pct = sum(w),
	      pos = min(position)) |>
    ungroup() |>
    group_by(readid) |>
    mutate(pct = pct/n_distinct(alig)) |>
    ungroup()

destinations_mate <- second_res |>
    mutate(readlen = map_int(positions, nrow)) |>
    unnest(cols = c(positions)) |>
    left_join(gene_annot, join_by(chr, between(position, start, end))) |>
    add_count(alig, readid, chr, position) |>
    mutate(w = (1/n)/readlen) |>
    group_by(alig, readid, gene_name) |>
    summarise(pct = sum(w),
	      pos = min(position)) |>
    ungroup() |>
    group_by(readid) |>
    mutate(pct = pct/n_distinct(alig)) |>
    ungroup()

destinations_df <- 
    bind_rows(first = destinations_read,
	      second = destinations_mate,
	      .id = "read") |>
    group_by(readid) |>
    mutate(pct = pct/sum(pct)) |>
    group_by(readid, gene_name) |>
    summarise(pct = sum(pct),
	      pos = median(pos)) |>
    ungroup()

# unmap?
unmap_df <- sam |>
    filter(flag %in% c(77, 141)) |>
    distinct(readid) |>
    mutate(gene_name = "unmap", pct = 1, pos = 1)

out <- bind_rows(destinations_df, unmap_df)

write_tsv(out, sprintf("./read_destinations/%s.tsv", sample_id))
