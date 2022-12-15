library(tidyverse)

phase1dir <- "/raid/genevol/nci_rnaseq/phase1/fastq" 
phase2dir <- "/raid/genevol/nci_rnaseq/phase2/fastq"

phase1 <- list.files(phase1dir, pattern = "*fastq.gz") |>
    tibble(fastq = _) |>
    mutate(sampleid = sub("^([^_]+).+$", "\\1", fastq)) |>
    distinct(sampleid)

phase2 <- list.files(phase2dir, pattern = "*fastq.gz") |>
    tibble(fastq = _) |>
    mutate(sampleid = sub("^([^_]+).+$", "\\1", fastq)) |>
    distinct(sampleid)

bind_rows("1" = phase1, "2" = phase2, .id = "phase") |>
    write_tsv("samples.txt", col_names=FALSE)


fastq_p1 <- list.files(phase1dir, pattern = "*fastq.gz") |>
    tibble(fastq = _) |>
    extract(fastq, c("sampleid", "lane", "r"), "(66K\\d+)_[^_]+_(L\\d+)_(R[12])_.*", remove = FALSE) |>
    group_by(sampleid, r) |>
    mutate(l = paste0("L", 1:n())) |>
    ungroup() |>
    unite("id", c(l, r), sep = "_") |>
    select(-lane) |>
    pivot_wider(names_from = id, values_from = fastq) |>
    mutate(reps = "biol rep 1")

fastq_p2 <- list.files(phase2dir, pattern = "*fastq.gz") |>
    tibble(fastq = _) |>
    extract(fastq, c("sampleid", "r"), "(66K\\d+)_merged_(R[12]).*", remove = FALSE) |>
    group_by(sampleid, r) |>
    mutate(l = paste0("L", 1:n())) |>
    ungroup() |>
    unite("id", c(l, r), sep = "_") |>
    pivot_wider(names_from = id, values_from = fastq) |>
    mutate(reps = "biol rep 2")

bind_rows(fastq_p1, fastq_p2) |>
    mutate(tissue = "PBMC",
	   title = paste(sampleid, tissue, reps, sep = ",")) |>
    select(sampleid, title, L1_R1, L1_R2, L2_R1, L2_R2) |>
    write_tsv("/raid/genevol/nci_rnaseq/fastq_list.tsv")

bind_rows(fastq_p1, fastq_p2) |>
    mutate(tissue = "PBMC",
	   title = paste(sampleid, tissue, reps, sep = ",")) |>
    select(sampleid, title, L1_R1, L1_R2, L2_R1, L2_R2) |>
    pivot_longer(-(sampleid:title), names_to = "i", values_to = "fastq") |>
    separate(i, c("lane", "r"), sep = "_") |>
    pivot_wider(names_from = r, values_from = fastq) |>
    drop_na() |>
    write_tsv("/raid/genevol/nci_rnaseq/fastq_list_long.tsv")

check_1 <- "/raid/genevol/nci_rnaseq/phase1/checksums.txt" |> 
    read_delim(delim = " ", col_names = c("checksum", "dummy", "fastq")) 

check_2 <- "/raid/genevol/nci_rnaseq/phase2/checksums.txt" |> 
    read_delim(delim = " ", col_names = c("checksum", "dummy", "fastq")) 

bind_rows("1" = check_1, "2" = check_2, .id = "batch") |>
    mutate(fastq = sub("fastq/", "", fastq)) |>
    arrange(batch, fastq) |>
    select(fastq, checksum) |>
    write_tsv("/raid/genevol/nci_rnaseq/checksums.tsv")

