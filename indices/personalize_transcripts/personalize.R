library(Biostrings)
library(hlaseqlib)
library(tidyverse)

# functions
replace_dots <- function(imgtseq) {
    
    m <- str_locate_all(imgtseq, "\\*\\.+\\*")[[1]]
    
    if (nrow(m) == 0) return(imgtseq)

    for (i in 1:nrow(m)) {
        st <- m[i, 1] + 1L
        en <- m[i, 2] - 1L
        substring(imgtseq, st, en) <- gsub("\\.", "*", substring(imgtseq, st, en))
    }
    
    return(imgtseq)
}

read_hla <- function(hla_locus) {
    
    hla_read_alignment(hla_locus, IMGTDB, imgtfile = "gen") %>%
	mutate(cds = map_chr(cds, replace_dots),
	       cds = gsub("\\*", "+", cds),
	       cds = gsub("\\.", "-", cds))
}


# HLA loci
loci <- c("A", "B", "C")

IMGTDB <- "/home/vitor/IMGTHLA"

# Gene annotations
annots <- "/home/vitor/gencode/gencode.v37.primary_assembly.annotation.gtf" %>% 
    read_tsv(comment = "#", col_names = FALSE, col_types = "c-cii-c-c")

gene_annots <- annots %>%
    filter(X3 == "gene") %>%
    transmute(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
              gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
              strand = X7) %>%
    filter(gene_name %in% paste0("HLA-", loci))

# match to genome
hlagen <- map_df(loci, read_hla) %>%
    filter(grepl("^\\+*[AGCT-]+\\+*$", cds)) %>%
    mutate(gene_name = sub("^([^*]+).+$", "HLA-\\1", allele)) %>%
    left_join(select(gene_annots, gene_name, strand), by = "gene_name") %>%
    select(gene_name, strand, allele, cds)

hlagen$cds[hlagen$strand == "-"] <- 
    as.character(reverseComplement(DNAStringSet(hlagen$cds[hlagen$strand == "-"])))

hlaset <- gsub("[+-]", "", hlagen$cds) %>%
    DNAStringSet() %>%
    setNames(hlagen$allele)

genome <- readDNAStringSet("/home/vitor/gencode/GRCh38.primary_assembly.genome.fa")

mhc <- subseq(genome["chr6 6"], 29e6, 32e6)

genomematch <- matchPDict(hlaset, mhc[[1]]) %>%
    as.list() %>%
    map_df(as.data.frame, .id = "allele") %>%
    as_tibble() %>%
    mutate_at(vars(start, end), ~`+`(., 29e6L - 1L)) %>%
    mutate(locus = sub("^([^*]+).+$", "\\1", allele)) %>%
    select(locus, allele, start, end, width)

write_tsv(genomematch, "genome_match.tsv")

# attribute positions to all alleles
hlagen_ref <- inner_join(hlagen, genomematch, by = "allele") %>%
    mutate(cdslen = map_int(cds, nchar))

hlagen_cds_split <- hlagen %>%
    mutate(i = map(cds, ~1:nchar(.)),
           cds = str_split(cds, "")) %>%
    select(gene_name, allele, i, cds) %>%
    unnest(c(i, cds))

hlagen_ref_cds_split <- hlagen_ref %>%
    mutate(i = map(cds, ~1:nchar(.)),
           cds = str_split(cds, "")) %>%
    select(gene_name, allele, i, cds) %>%
    unnest(c(i, cds))

sequence_status <- hlagen_ref_cds_split %>%
    mutate(status = ifelse(cds == "+", "unknown", "known")) %>%
    select(gene_name, i, status)

hla_ref_pos <- hlagen_ref %>%
    mutate(pos = map2(start, end, ~.x:.y),
           ix = map(lengths(pos), seq_len)) %>%
    select(gene_name, allele, ix, pos) %>%
    unnest(c(ix, pos))

hla_map <- hlagen_ref_cds_split %>%
    filter(cds != "+", cds != "-") %>%
    group_by(gene_name, allele) %>%
    mutate(ix = 1:n()) %>%
    ungroup() %>%
    left_join(hla_ref_pos, by = c("gene_name", "allele", "ix")) %>%
    select(gene_name, i, pos) %>%
    group_by(gene_name) %>%
    nest() %>%
    ungroup() %>%
    left_join(select(hlagen_ref, gene_name, cdslen), by = "gene_name") %>%
    mutate(data = map2(data, cdslen, ~complete(.x, i = 1:.y))) %>%
    select(-cdslen) %>%
    unnest(data) %>%
    left_join(sequence_status, by = c("gene_name", "i")) %>%
    group_by(gene_name) %>%
    fill(pos) %>%
    ungroup() %>%
    mutate(pos = ifelse(status == "unknown", NA, pos))

positions_all <- hlagen_cds_split %>%
    right_join(hla_map, by = c("gene_name", "i")) %>%
    filter(cds != "+", status  == "known") %>%
    select(-status)


# transcripts
hla_transcripts_annot <- annots %>%
    filter(X3 == "transcript") %>%
    transmute(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	      gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"), 
	      tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+")) %>%
    filter(gene_name %in% paste0("HLA-", c("A", "B", "C")))

exon_annots <- annots %>%
    filter(X1 == "chr6", X3 == "exon") %>%
    transmute(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
              tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
	      tx_type = str_extract(X9, "(?<=transcript_type\\s\")[^\"]+"),
	      start = X4, 
              end = X5) %>%
    inner_join(gene_annots, by = "gene_id") %>%
    select(gene_id, gene_name, tx_id, strand, start, end) %>%
    mutate(sq = map2_chr(start, end, ~as.character(subseq(genome["chr6 6"], .x, .y)))) %>%
    arrange(gene_name, tx_id, start)

exon_positions <- exon_annots %>%
    mutate(pos = map2(start, end, ~.x:.y),
	   sq = str_split(sq, "")) %>%
    select(gene_name, tx_id, sq, pos) %>%
    unnest(c(pos, sq))


# HLA genotypes
hla_genotypes <- 
    read_tsv("../../analysis/genos.tsv") %>%
    select(-ends_with("_ok")) %>%
    pivot_longer(-(1:3), names_to = "method", values_to = "allele") %>%
    mutate(method = sub("^allele_", "", method),
	   allele = sub("G$", "", allele)) %>%
    filter(!is.na(allele)) %>%
    separate_rows(allele, sep = ";") %>%
    mutate(allele_full = map_chr(allele, 
				 ~first(grep(., hlagen$allele, fixed = TRUE, value = TRUE)))) %>%
    group_by(sampleid, gene_name, i) %>%
    filter( all(is.na(allele_full)) | !is.na(allele_full)  ) %>%
    ungroup() %>%
    distinct(sampleid, gene_name, i, allele = allele_full)

write_tsv(hla_genotypes, "./hla_genotypes_index.tsv")


hla_alleles <- hla_genotypes %>%
    distinct(gene_name, allele = allele_full) %>%
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

results$sq[results$strand == "-"] <- 
    as.character(reverseComplement(DNAStringSet(results$sq[results$strand == "-"])))

hla_index <- DNAStringSet(results$sq) %>%
    setNames(paste(results$tx_id, results$allele, sep = "_"))

writeXStringSet(hla_index, "personalized_transcripts.fa")
