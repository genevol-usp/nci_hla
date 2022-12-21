library(tidyverse)

transcript_annot <- read_tsv("./transcript_annot.tsv")

index <- Biostrings::readDNAStringSet("../make_simulation_data/simulation_index.fa")

simul_counts <- read_tsv("../make_simulation_data/phenotypes.tsv") |>
    add_column(txid = names(index), .before = 1) |>
    mutate(txid = sub("^([^_]+).*$", "\\1", txid)) |>
    left_join(transcript_annot, by = c("txid" = "transcript_id")) |>
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) |>
    select(gene_id, gene_name, txid, starts_with("66K")) |>
    pivot_longer(starts_with("66K"), names_to = "sampleid") |>
    group_by(sampleid, gene_name) |>
    summarise(value = sum(value)) |>
    ungroup()

read_origin <- read_tsv("read_origins_salmon.tsv")

samples <- unique(read_origin$sampleid)

read_dest <- "./read_destinations/%s_salmon.tsv" |>
    sprintf(samples) |>
    setNames(samples) |>
    map_df(read_tsv, .id = "sampleid")

tdf <- 
    inner_join(read_dest, read_origin, 
	       by = c("sampleid", "readid"), 
	       suffix = c("_dest", "_origin")) |> 
    select(sampleid, readid, origin = gene_name_origin, pos_orig = origin_read,
	   dest = gene_name_dest, pos_dest = pos, pct) |>
    group_by(sampleid, origin, dest) |>
    summarise(pct = sum(pct)) |>
    group_by(sampleid, dest) |>
    summarise(totalmap = sum(pct)) |>
    ungroup() |>
    filter(origin %in% c("HLA-A", "HLA-B", "HLA-C") | 
	   dest %in% c("HLA-A", "HLA-B", "HLA-C")) |>
    left_join(simul_counts, by = c("sampleid", "origin" = "gene_name"))

# add counts for reads originated elsewhere than HLA-A, -B, -C
# interpret these counts as "gain" in respect to simulated
tmp <- tdf |>
    filter(is.na(value)) |>
    select(-value) |>
    left_join(simul_counts, by = c("sampleid", "dest" = "gene_name"))

final <- filter(tdf, !is.na(value)) |>
    bind_rows(tmp) |>
    mutate(w_by_simul = pct/value,
	   w_by_mapped = pct/totalmap) |>
    select(sampleid, origin, dest, w_by_simul, w_by_mapped) |> 
    filter(!is.na(dest)) # filter out soft clips

write_tsv(final, "./origin_dest_results_salmon.tsv")
