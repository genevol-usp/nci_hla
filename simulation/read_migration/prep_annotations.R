library(tidyverse)

# Gencode annotations
annotations <- 
    "/raid/genevol/gencode/gencode.v37.primary_assembly.annotation.gtf" |> 
    read_tsv(comment = "#", col_types = "c-cii-c-c",
	     col_names = c("chr", "feature", "start", "end", "strand", "info"))

gene_annot <- annotations |>
    filter(feature == "gene") |>
    mutate(gene_id = str_extract(info, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+")) |>
    select(-feature, -info, -strand)

write_tsv(gene_annot, "gene_coords.tsv")

#exon_annot <- annotations |>
#    filter(feature == "exon") |>
#    mutate(transcript_id = str_extract(info, "(?<=transcript_id\\s\")[^\"]+"),
#	   gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+")) |>
#    select(-feature, -info) 
#

transcript_annot <- annotations |>
    filter(feature == "transcript") |>
    mutate(transcript_id = str_extract(info, "(?<=transcript_id\\s\")[^\"]+"),
	   gene_id = str_extract(info, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+")) |>
    select(transcript_id, gene_id, gene_name, strand)
    
write_tsv(transcript_annot, "./transcript_annot.tsv")
