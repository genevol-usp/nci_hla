library(tidyverse)

annots <- "~/gencode/gencode.v37.primary_assembly.annotation.gtf" %>% 
    read_tsv(comment = "#", col_names = FALSE, col_types = "c-cii-c-c")

gene_annots <- annots %>%
    filter(X3 == "gene") %>%
    transmute(chr = X1,
	      gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	      gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	      start = X4, end = X5, strand = X7)

tx_annots <- annots %>%
    filter(X3 == "transcript") %>%
    transmute(chr = X1,
	      gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	      gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	      tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
	      tx_type = str_extract(X9, "(?<=transcript_type\\s\")[^\"]+"),  
	      start = X4, end = X5, strand = X7)

hla_exon_annots <- annots %>%
    filter(X3 == "exon") %>%
    transmute(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	      gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	      tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
	      exon_id = str_extract(X9, "(?<=exon_id\\s\")[^\"]+"),
	      start = X4, end = X5, strand = X7) %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C"))

write_tsv(gene_annots, "./gene_annotation_df.tsv")
write_tsv(tx_annots, "./transcript_annotation_df.tsv")
write_tsv(hla_exon_annots, "./hla_exon_annots.tsv")

