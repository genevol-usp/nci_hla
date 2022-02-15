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
