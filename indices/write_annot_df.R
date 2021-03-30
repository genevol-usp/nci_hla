library(tidyverse)

annots <- "~/gencode/gencode.v37.primary_assembly.annotation.gtf" %>% 
    read_tsv(comment = "#", col_names = FALSE, col_types = "c-cii-c-c") %>%
    filter(X3 == "transcript") %>%
    transmute(chr = X1,
	      gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	      gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	      tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
	      start = X4, end = X5, strand = X7)

write_tsv(annots, "./transcript_annotation_df.tsv")

