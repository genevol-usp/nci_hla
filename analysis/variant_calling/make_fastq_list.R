library(tidyverse)

phase1 <- "/raid/genevol/nci_rnaseq/phase1/fastq" %>%
    list.files(full.names = TRUE)

phase2 <- "/raid/genevol/nci_rnaseq/phase2/fastq" %>%
    list.files(full.names = TRUE)

fastq_df <- 
    bind_rows("1" = tibble(path = phase1),
	      "2" = tibble(path = phase2),
	      .id = "phase") %>%
    mutate(sampleid = str_remove(basename(path), "_R[12].+$"),
	   strand = str_extract(path, "(R[12])")) %>%
    pivot_wider(names_from = "strand", values_from = path) %>%
    select(-phase)

write_tsv(fastq_df, "./fastq_list.txt", col_names = FALSE)

bam_df <- fastq_df %>%
    mutate(subject_id = sub("^([^_]+).+$", "\\1", sampleid),
	   bam = paste0("/scratch/vitor/results/gatk/", sampleid, "_recal.bam")) %>%
    group_by(subject_id) %>%
    summarise(bam = paste(bam, collapse = ",")) %>%
    ungroup()

write_tsv(bam_df, "./bam_list.txt", col_names = FALSE)

tibble(fastq = phase1) %>%
    mutate(sampleid = str_remove(basename(fastq), "_R[12].+$"),
	   subject_id = sub("^(66K\\d+)_.+$", "\\1", sampleid)) %>%
    distinct(subject_id, sampleid) %>%
    mutate(path = paste0("/scratch/vitor/results/star/pass_wasp/", sampleid, ".aseinput.bam")) %>%
    group_by(subject_id) %>%
    summarise(path = paste(path, collapse = ",")) %>%
    ungroup() %>%
    write_tsv("./bam_ase_list.txt", col_names = FALSE)
