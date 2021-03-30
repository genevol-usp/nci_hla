library(Biostrings)
library(tidyverse)

annots <- "~/gencode/gencode.v37.primary_assembly.annotation.gtf" %>% 
    read_tsv(comment = "#", col_names = FALSE, col_types = "c-cii-c-c") %>%
    filter(X3 == "transcript") %>%
    transmute(chr = X1,
	      gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	      tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
	      start = X4, end = X5)

#transcripts <- readDNAStringSet("../make_simulation_data/simulation_index.fa") %>%
#    .[grep("_[ABC]\\*", names(.))]
#
#ids <- sub("^([^_]+).*$", "\\1", names(transcripts)) %>%
#    unique()

bed <- annots %>% 
    #filter(tx_id %in% ids) %>%
    filter(chr == "chr6",
	   gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    group_by(chr, gene_name) %>%
    summarise(start = min(start),
	      end = max(end)) %>%
    ungroup() %>%
    select(chr, start, end, gene_name)

write_tsv(bed, "hla.bed", col_names = FALSE)
