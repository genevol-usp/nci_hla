library(Biostrings)
library(tidyverse)

loci <- c("A", "B", "C")

annots <- "/raid/genevol/gencode/gencode.v37.primary_assembly.annotation.gtf" %>% 
    read_tsv(comment = "#", col_names = FALSE, col_types = "c-cii-c-c")

gene_annots <- annots %>%
    filter(X3 == "gene") %>%
    transmute(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
              gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
              strand = X7) %>%
    filter(gene_name == "HLA-A")

genome <- readDNAStringSet("/raid/genevol/gencode/GRCh38.primary_assembly.genome.fa")

positions_all <- read_rds("./hla_allele_genome_map.rds") %>%
    filter(gene_name == "HLA-A")

# check last position
positions_all %>%
    filter(allele == "A*01:01:01:01") %>%
    tail()

# transcripts
hla_transcripts_annot <- annots %>%
    filter(X3 == "transcript") %>%
    transmute(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	      gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"), 
	      tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+")) %>%
    filter(gene_name == "HLA-A")

exon_annots <- annots %>%
    filter(X1 == "chr6", X3 == "exon") %>%
    transmute(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
              tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
	      tx_type = str_extract(X9, "(?<=transcript_type\\s\")[^\"]+"),
	      start = X4, 
              end = X5) %>%
    inner_join(gene_annots, by = "gene_id") %>%
    select(gene_id, gene_name, tx_id, strand, start, end) %>%
    mutate(end = ifelse(end > 29945765, 29945765, end)) %>%
    mutate(sq = map2_chr(start, end, ~as.character(subseq(genome["chr6 6"], .x, .y)))) %>%
    arrange(gene_name, tx_id, start)

exon_positions <- exon_annots %>%
    mutate(pos = map2(start, end, ~.x:.y),
	   sq = str_split(sq, "")) %>%
    select(gene_name, tx_id, sq, pos) %>%
    unnest(c(pos, sq))


# HLA genotypes
hla_genotypes <- read_tsv("./hla_genotypes_index.tsv") %>%
    filter(gene_name == "HLA-A")

hla_alleles <- hla_genotypes %>%
    distinct(gene_name, allele) %>%
    arrange(gene_name, allele)

# personalize transcripts
personalize_transc <- function(x_allele, y_transcript) {

    filter(positions_all, allele == x_allele) %>%
    left_join(filter(exon_positions, tx_id == y_transcript), ., 
	      by = c("gene_name", "pos")) %>%
    mutate(final_seq = ifelse(is.na(cds), sq, cds)) %>%
    group_by(tx_id) %>%
    summarise(final_seq = paste(final_seq, collapse = "")) %>%
    ungroup() %>%
    pull(final_seq) %>%
    gsub("-", "", .)
}

## pairs of transcripts and HLA alleles to run personalization on
apply_set <- inner_join(hla_alleles,
			distinct(exon_positions, gene_name, tx_id),
			by = "gene_name")

## run personalization
results <- apply_set %>%
    mutate(sq = map2_chr(allele, tx_id, personalize_transc)) %>%
    left_join(gene_annots, by = "gene_name") 

hla_index <- DNAStringSet(results$sq) %>%
    setNames(paste(results$tx_id, results$allele, sep = "_"))

hla_index_ori <- readDNAStringSet("./personalized_transcripts.fa")

all_seqs <- c(hla_index_ori, hla_index)

all_seqs_df <- tibble(id = names(all_seqs), s = as.character(all_seqs)) %>%
    distinct(id, s)

all_seqs_a <- all_seqs_df %>%
    filter(grepl("_A\\*", id)) %>%
    mutate(len = nchar(s)) %>%
    group_by(id) %>%
    arrange(id, len) %>%
    mutate(lenclass = case_when(n() == 1 ~ "u",
				n() == 2 & len == min(len) ~ "s",
				n() == 2 & len == max(len) ~ "l",
				TRUE ~ NA_character_)) %>%
    ungroup() %>%
    select(-len) %>%
    unite(id, c("id", "lenclass"), sep = ".")

all_seqs_final <- all_seqs_df %>%
    filter(!grepl("_A\\*", id)) %>%
    bind_rows(all_seqs_a, .)

out <- DNAStringSet(all_seqs_final$s) %>% 
    setNames(all_seqs_final$id)

writeXStringSet(out, "personalized_transcripts_plusShort.fa")
