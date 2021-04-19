library(Biostrings)
library(tidyverse)

# sample IDs
samples <- readLines("./samples.txt")

# compute sequence divergence of personalized transcripts
# in regard to original sequences
transcripts <- "./make_simulation_data/simulation_index.fa" %>%
    readDNAStringSet()

pers_isoforms <- transcripts[grepl("\\*", names(transcripts))] 

hla_transcripts <- names(pers_isoforms) %>%
    str_split("_") %>%
    map(1) %>%
    unlist() %>%
    unique()

ref_isoforms <- "../indices/gencode.transcripts.fa" %>%
    readDNAStringSet() %>%
    .[hla_transcripts]

iso_df <- tibble(pers = names(pers_isoforms)) %>%
    separate(pers, c("tx", "allele"), sep = "_", remove = FALSE)

submat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1)

calc_identity <- function(pers_id, ref_id) {
    x <- pers_isoforms[pers_id]
    y <- ref_isoforms[ref_id]
    
    pairwiseAlignment(x, y, substitutionMatrix = submat, gapOpening = -1, gapExtension = 0) %>%
        pid()
}

identity_df <- iso_df %>%
    mutate(ident = map2_dbl(pers, tx, calc_identity))

write_rds(identity_df, "./plot_data/simulated_transcript_identity.rds")


# transcript annotations

annots <- read_tsv("../indices/transcript_annotation_df.tsv") %>%
    select(gene_id, gene_name, tx_id)

hla_annots <- 
    tibble(pers = grep("\\*", names(transcripts), value = TRUE),
	   tx_id = sub("^([^_]+).*$", "\\1", pers)) %>%
    left_join(annots, by = "tx_id") %>%
    select(gene_id, gene_name, tx_id = pers)

all_annots <- bind_rows(annots, hla_annots)


# True expression
ground_truth <- "./make_simulation_data/phenotypes.tsv" %>%
    read_tsv() %>%
    add_column(tx_id = names(transcripts), .before = 1) %>%
    left_join(all_annots, by = "tx_id") %>%
    pivot_longer(starts_with("66K"), names_to = "sampleid", values_to = "true_counts") %>%
    select(sampleid, gene_id, gene_name, tx_id, true_counts)

## gene-level true expression
ground_truth_genes <- ground_truth %>%
    group_by(sampleid, gene_id, gene_name) %>%
    summarise(true_counts = sum(true_counts)) %>%
    ungroup()

# Compute weighted distances
dist_df <- ground_truth %>%
    filter(grepl("\\*", tx_id)) %>%
    filter(true_counts > 0) %>%
    arrange(sampleid) %>%
    select(sampleid, gene_id, gene_name, tx_id, true_counts) %>%
    left_join(select(identity_df, tx_id = pers, ident), by = "tx_id") %>%
    group_by(sampleid, gene_name) %>%
    summarise(ident = weighted.mean(ident, w = true_counts)) %>%
    ungroup() %>%
    mutate(wdist = (100L - ident)/100L) %>%
    select(sampleid, gene_name, wdist)

write_rds(dist_df, "./plot_data/weighted_distances.rds")

# Quantifications
## Salmon Ref
salmon_ref <- sprintf("./salmon/quant/%s/quant.sf", samples) %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "sampleid") %>%
    left_join(annots, by = c("Name" = "tx_id")) %>%
    group_by(sampleid, gene_id, gene_name) %>%
    summarise(counts = sum(NumReads)) %>%
    ungroup() %>%
    full_join(ground_truth_genes, by = c("sampleid", "gene_id", "gene_name")) %>%
    mutate_at(vars(counts, true_counts), ~replace_na(., 0))

## Salmon pers
salmon_pers <- sprintf("./salmon-pers/quant/%s/quant.sf", samples) %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "sampleid") %>%
    left_join(all_annots, by = c("Name" = "tx_id")) %>%
    group_by(sampleid, gene_id, gene_name) %>%
    summarise(counts = sum(NumReads)) %>%
    ungroup() %>%
    full_join(ground_truth_genes, by = c("sampleid", "gene_id", "gene_name")) %>%
    mutate_at(vars(counts, true_counts), ~replace_na(., 0))

## HLApers
hlapers_annot <- sprintf("./hlapers/genotypes/%s_genotypes.tsv", samples) %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "sampleid") %>%
    mutate(allele = sub("^([^-]+).*$", "\\1", allele),
	   gene_name = paste0("HLA-", locus)) %>%
    left_join(distinct(annots, gene_id, gene_name), by = "gene_name") %>%
    distinct(gene_id, gene_name, tx_id = allele)

hlapers <- sprintf("./hlapers/quant/%s_quant/quant.sf", samples) %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "sampleid") %>%
    mutate(Name = sub("^([^-]+).*$", "\\1", Name)) %>%
    left_join(bind_rows(annots, hlapers_annot), by = c("Name" = "tx_id")) %>%
    group_by(sampleid, gene_id, gene_name) %>%
    summarise(counts = sum(NumReads)) %>%
    ungroup() %>%
    full_join(ground_truth_genes, by = c("sampleid", "gene_id", "gene_name")) %>%
    mutate_at(vars(counts, true_counts), ~replace_na(., 0))

## hla-mapper::rna
hlamapper <- sprintf("./hla-mapper/quant/%s_counts.txt", samples) %>%
    setNames(samples) %>%
    map_df(~read_tsv(., comment = "#", col_types = "c----dd") %>% 
	   select(gene_id = 1, len = 2, counts = 3), 
           .id = "sampleid") %>%
    left_join(distinct(annots, gene_id, gene_name), by = "gene_id") %>%
    full_join(ground_truth_genes, by = c("sampleid", "gene_id", "gene_name")) %>%
    mutate_at(vars(counts:true_counts), ~replace_na(., 0)) %>%
    select(sampleid, gene_id, gene_name, counts, true_counts)

quant_df <- 
    bind_rows("Salmon ref genome" = salmon_ref,
	      "Salmon personalized" = salmon_pers,
	      "HLApers" = hlapers,
	      "hla-mapper::rna" = hlamapper,
	      .id = "method") %>%
    mutate(method = factor(method, levels = unique(method)))


# Assess accuracy

## Compute TPM/TRUE for HLA
count_rates <- quant_df %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    mutate(rate_counts = counts/true_counts) %>%
    select(method, sampleid, gene_id, gene_name, rate_counts) %>%
    arrange(method, sampleid, gene_name) %>%
    mutate(gene_name = factor(gene_name, levels = unique(gene_name)))

write_rds(count_rates, "./plot_data/count_rates.rds")

# hla-mapper coverage
#
### exon coordinates
#hlamatch <- read_tsv("./make_simulation_data/genome_match.tsv") %>%
#    mutate(pos = case_when(locus == "A" ~ map2(start, end, ~.x:.y),
#			   locus == "B" | locus == "C" ~ map2(start, end, ~.y:.x))) %>%
#    select(locus, allele, pos) %>%
#    unnest(pos) %>%
#    group_by(locus, allele) %>%
#    mutate(i = seq_len(n())) %>%
#    ungroup()
#
#hla_coords <- hlamatch$locus %>%
#    map_df(~hlaseqlib::hla_read_alignment(., imgtdb = "~/IMGTHLA", imgtfile = "gen", by_exon = TRUE) %>%
#	   filter(allele %in% hlamatch$allele) %>%
#	   mutate(cds = gsub("\\.|\\*", "", cds),
#		  cds = str_split(cds, "")) %>%
#	   unnest(cds)) %>%
#    group_by(allele) %>%
#    mutate(i = seq_len(n())) %>%
#    ungroup()
#
#hla_exon_coords <- left_join(hla_coords, hlamatch, by = c("allele", "i")) %>%
#    mutate(gene_name = paste0("HLA-", locus)) %>%
#    unite(feature, c("feature", "idx_grp"), sep = " ") %>%
#    select(gene_name, feature, i, pos)
#
#bed <- "./hla-mapper/hla.bed" %>%
#    read_tsv(col_names = FALSE) %>%
#    select(gene_name = X4, start = X2, end = X3) %>%
#    mutate(pos = map2(start, end, ~.x:.y)) %>%
#    select(-start, -end) %>%
#    unnest(pos)
#
#coverage_df <- sprintf("./hla-mapper/coverage/%s.cov", samples) %>%
#    setNames(samples) %>%
#    map_df(~read_tsv(., col_names = FALSE), .id = "sampleid") %>%
#    select(sampleid, pos = X2, cov = X3) %>%
#    left_join(bed, by = "pos") %>%
#    inner_join(hla_exon_coords, by = c("gene_name", "pos"))
#
#write_rds(coverage_df, "./plot_data/coverage.rds")
