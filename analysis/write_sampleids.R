library(tidyverse)

phase1dir <- "/raid/genevol/nci_rnaseq/phase1/fastq" 
phase2dir <- "/raid/genevol/nci_rnaseq/phase2/fastq"

phase1 <- list.files(phase1dir) %>%
    tibble(fastq = .) %>%
    mutate(sampleid = sub("^([^_]+).+$", "\\1", fastq)) %>%
    distinct(sampleid)

phase2 <- list.files(phase2dir) %>%
    tibble(fastq = .) %>%
    mutate(sampleid = sub("^([^_]+).+$", "\\1", fastq)) %>%
    distinct(sampleid)

bind_rows("1" = phase1, "2" = phase2, .id = "phase") %>%
    write_tsv("samples.txt", col_names=FALSE)
