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
    select(gene_id, gene_name, tx_id, tx_type, strand, start, end) %>%
    filter(tx_type == "protein_coding") %>%
    select(-tx_type) %>%
    mutate(sq = map2_chr(start, end, ~as.character(subseq(genome["chr6 6"], .x, .y)))) %>%
    arrange(gene_name, tx_id, start)

hla_transcripts <- distinct(exon_annots, gene_name, tx_id)

## select transcripts which explain 90% of the expression in the real data
nci_samples <- read_tsv("../../analysis/samples.txt", col_names = c("tp", "sampleid")) %>%
    filter(tp == 1)

nci_expression <- "../../analysis/salmon/quant/%s_t1/quant.sf" %>%
    sprintf(nci_samples$sampleid) %>%
    setNames(nci_samples$sampleid) %>%
    map_df(read_tsv, .id = "sampleid") %>%
    group_by(sampleid) %>%
    mutate(scaled_counts = as.integer(round(NumReads/sum(NumReads) * 3e7))) %>%
    ungroup()

hla_expression <- inner_join(nci_expression, hla_transcripts, by = c("Name" = "tx_id")) 

hla_selected_transcripts <- hla_expression %>%
    group_by(sampleid, gene_name) %>%
    mutate(avg_tpm = (TPM/sum(TPM)) * 100) %>%
    group_by(gene_name, tx_id = Name) %>%
    summarise(avg_tpm = mean(avg_tpm)) %>%
    ungroup() %>%
    arrange(gene_name, -avg_tpm) %>%
    group_by(gene_name) %>%
    mutate(cumexp = cumsum(avg_tpm),
	   i = cumsum(cumexp >= 90)) %>%
    ungroup() %>%
    filter(i <= 1L)

exon_positions <- exon_annots %>%
    filter(tx_id %in% unique(hla_selected_transcripts$tx_id)) %>%
    mutate(pos = map2(start, end, ~.x:.y),
	   sq = str_split(sq, "")) %>%
    select(gene_name, tx_id, sq, pos) %>%
    unnest(c(pos, sq))


# HLA genotypes
hla_genotypes <- 
    "/raid/genevol/nci_rnaseq/phase1/HLA-A, -B, -C expression levels for Pat.xlsx" %>%
    readxl::read_xlsx() %>%
    select(sampleid = 3, starts_with("Class")) %>%
    pivot_longer(-1, names_to = "hap", values_to = "allele") %>%
    mutate(hap = sub("Class 1", "", hap)) %>%
    separate(hap, c("gene_name", "i"), sep = " ") %>%
    mutate(i = parse_number(i)) %>%
    extract(allele, c("f1", "f2", "f3"), "(\\d{2})(\\d{2})(\\d*)") %>%
    mutate(f3 = if_else(f3 != "", str_pad(f3, 2, pad = "0"), f3)) %>%
    unite(allele, c("f1", "f2", "f3"), sep = ":") %>%
    mutate(allele = sub(":$", "", allele),
	   allele = paste(gene_name, allele, sep = "*"),
	   allele = ifelse(grepl("NA", allele), NA, allele),
	   allele = recode(allele, "C*17:00" = "C*17"),
	   allele = ifelse(sampleid == "66K00241" & allele == "A*03:01:01", "A*30:01:01", allele),
	   gene_name = paste0("HLA-", gene_name),
	   inferred = map_chr(allele, ~first(grep(., hlagen$allele, fixed = TRUE, value = TRUE)))) %>%
    select(sampleid, gene_name, i, allele = inferred) %>%
    group_by(sampleid) %>%
    filter(!any(is.na(allele))) %>%
    ungroup()

## choose 50 individuals
set.seed(10)
simul_df <- hla_genotypes %>%
    filter(sampleid %in% hla_expression$sampleid) %>%
    group_by(sampleid) %>%
    nest() %>%
    ungroup() %>%
    arrange(sampleid) %>%
    sample_n(50) %>%
    unnest(c(data))

write_tsv(simul_df, "./hlagenotypes.tsv")

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
apply_set <- inner_join(distinct(simul_df, gene_name, allele),
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

# make count matrix
genome_background_quants <- nci_expression %>%
    filter(sampleid == "66K00003") %>%
    select(tx_id = Name, readcount = scaled_counts) %>%
    filter(! tx_id %in% hla_transcripts_annot$tx_id)

hla_ground_truth <- hla_expression %>%
    filter(Name %in% unique(hla_selected_transcripts$tx_id)) %>%
    select(sampleid, gene_name, tx_id = Name, readcount = scaled_counts) %>%
    arrange(sampleid, gene_name, tx_id) %>%
    inner_join(simul_df, by = c("sampleid", "gene_name")) %>%
    group_by(sampleid, tx_id) %>%
    mutate(readcount = readcount/2L) %>%
    ungroup() %>%
    unite(tx_id, c("tx_id", "allele"), sep = "_") %>%
    group_by(sampleid, tx_id) %>%
    summarise(readcount = sum(readcount)) %>%
    ungroup()

final_quants <- hla_ground_truth %>%
    split(.$sampleid) %>%
    map(~select(., -sampleid)) %>%
    map_df(~bind_rows(genome_background_quants, .), .id = "sampleid")

pheno_matrix <- final_quants %>%
    pivot_wider(names_from = sampleid, values_from = readcount) %>%
    mutate_at(vars(-tx_id), ~replace_na(., 0))

# index for simulation
index_transcriptome <- "../../indices/gencode.transcripts.fa" %>%
    readDNAStringSet()

index <- c(index_transcriptome, hla_index)[pheno_matrix$tx_id]

# save output
write_tsv(select(pheno_matrix, -1), "phenotypes.tsv")
writeXStringSet(index, "simulation_index.fa")

# compute mean and sd of fragment length distribution
fld <- scan("../../analysis/salmon/quant/66K00003_t1/libParams/flenDist.txt")
    setnames(sprintf(
    map(

fld_params <- tibble(mu = sum(1:length(fld) * fld)) %>%
    mutate(sigma = sqrt(sum( ( (1:length(fld) - mu)^2 ) * (fld) )))

write_tsv(fld_params, "fld_params.tsv")
