library(tidyverse)

annot <- 
    file.path("/raid/genevol/gencode/gencode.v37.primary_assembly.annotation.gtf") %>%
    read_tsv(comment = "#", col_types = "c-cii-c-c",
	     col_names = c("chr", "feature", "start", "end", "strand", "info")) %>%
    mutate(chr = factor(chr, levels = str_sort(unique(chr), numeric = TRUE)))

exon_hla <- annot %>%
    filter(feature == "exon") %>%
    mutate(gene_id = str_extract(info, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+"),
	   transcript_id = str_extract(info, "(?<=transcript_id\\s\")[^\"]+")) %>%
    select(chr, start, end, strand, gene_id, gene_name, transcript_id) %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C", "B2M"))

intron_annot <- exon_hla %>%
    group_by(gene_id) %>%
    mutate(start_gene = min(start), end_gene = max(end)) %>%
    ungroup() %>%
    arrange(chr, start_gene, start) %>%
    group_by(gene_id, gene_name, transcript_id, start_gene, end_gene) %>%
    summarise(data = bind_cols(tibble(start = c(unique(start_gene), end)),
			       tibble(end = c(start, unique(end_gene))))) %>%
    ungroup() %>%
    unnest(cols = c(data)) %>%
    mutate(start_i = ifelse(start == end | start == start_gene, start, start + 1L),
	   end_i = ifelse(start == end | end == end_gene, end, end - 1L)) %>%
    filter(start != end) %>%
    select(gene_id, gene_name, transcript_id, start = start_i, end = end_i)

annot_out <- 
    bind_rows("exon" = select(exon_hla, gene_id, gene_name, transcript_id, start, end),
              "intron" = intron_annot,
              .id = "feature") %>%
    arrange(gene_name, gene_id, transcript_id, start)

write_rds(annot_out, "./plot_data/transcripts_annot.rds")

gene_hla <- annot %>%
    filter(feature == "gene") %>%
    filter(grepl("HLA-[ABC]", info)) %>%
    mutate(gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+")) %>%
    select(chr, gene_name, start, end)


gene_hla %>%
    select(chr, start, end) %>%
    mutate(start = start - 1L) %>%
    write_tsv("./plot_data/hla.bed", col_names=FALSE)
