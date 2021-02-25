library(Biostrings)
library("BSgenome.Hsapiens.NCBI.GRCh38")
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

# Gene annotations -> strand
annots <- "/home/vitor/hisat2/grch38_snp_tran/Homo_sapiens.GRCh38.99.gtf" %>% 
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

genomematch <- matchPDict(hlaset, Hsapiens$"6"[29e6:32e6]) %>%
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
exon_annots <- annots %>%
    filter(X1 == 6) %>%
    filter(X3 == "exon") %>%
    transmute(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
              tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
	      tx_biotype = str_extract(X9, "(?<=transcript_biotype\\s\")[^\"]+"),
	      start = X4, 
              end = X5) %>%
    inner_join(gene_annots, by = "gene_id") %>%
    select(gene_id, gene_name, tx_id, tx_biotype, strand, start, end) %>%
    filter(tx_biotype == "protein_coding") %>%
    filter(tx_id != "ENST00000640615") %>% #B isoform with exon from C; computationally predicted
    select(-tx_biotype) %>%
    mutate(sq = map2_chr(start, end, ~substring(as.character(Hsapiens$"6"), .x, .y))) %>%
    arrange(gene_name, tx_id, start)

hla_transcripts <- distinct(exon_annots, gene_name, tx_id)

## select transcripts which explain 90% of the expression in the real data
nci_samples <- read_tsv("../../analysis/samples.txt", col_names = c("tp", "sampleid")) %>%
    filter(tp == 1)

hla_expression <- "../../analysis/salmon/quant/%s_t1/quant.sf" %>%
    sprintf(nci_samples$sampleid) %>%
    setNames(nci_samples$sampleid) %>%
    map_df(~read_tsv(.) %>% 
	   inner_join(hla_transcripts, by = c("Name" = "tx_id")), 
           .id = "sampleid")

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
    filter(sampleid %in% sample(unique(sampleid), 50))

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
genome_background_quants <- "../../analysis/salmon/quant/66K00003_t1/quant.sf" %>%
    read_tsv() %>%
    select(tx_id = Name, readcount = NumReads) %>%
    filter(! tx_id %in% hla_transcripts$tx_id)

hla_ground_truth <- hla_expression %>%
    filter(Name %in% unique(hla_selected_transcripts$tx_id)) %>%
    select(sampleid, gene_name, tx_id = Name, readcount = NumReads) %>%
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
    map_df(~bind_rows(genome_background_quants, .) %>%
	   fill(sampleid, .direction = "up")) %>%
    group_by(sampleid) %>%
    mutate(readcount = as.integer(round(readcount/sum(readcount) * 3e7))) %>%
    ungroup()

pheno_matrix <- final_quants %>%
    pivot_wider(names_from = sampleid, values_from = readcount) %>%
    mutate_at(vars(-tx_id), ~replace_na(., 0))

# index for simulation
index_transcriptome <- "/raid/genevol/geuvadis/salmon/ensembl.transcripts.fa" %>%
    readDNAStringSet()

index <- c(index_transcriptome, hla_index)[pheno_matrix$tx_id]

# save output
write_tsv(select(pheno_matrix, -1), "phenotypes.tsv")
writeXStringSet(index, "simulation_index.fa")

# compute mean and sd of fragment length distribution
fld <- scan("../../analysis/salmon/quant/66K00003_t1/libParams/flenDist.txt")

fld_params <- tibble(mu = sum(1:length(fld) * fld)) %>%
    mutate(sigma = sqrt(sum( ( (1:length(fld) - mu)^2 ) * (fld) )))

write_tsv(fld_params, "fld_params.tsv")