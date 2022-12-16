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

transc_lens <- 
    tibble(tx_id = names(transcripts),
           len = width(transcripts)) %>%
    mutate(eff_len = len - 261 + 1,
           eff_len = ifelse(eff_len < 0, 0, eff_len)) %>%
    select(tx_id, eff_len)

ground_truth_tx_tpms <- ground_truth %>%
    filter(true_counts > 0) %>%
    left_join(transc_lens) %>%
    group_by(sampleid) %>%
    mutate(rate = log(true_counts) - log(eff_len),
           rate = ifelse(eff_len == 0, 0, rate),
           denom = log(sum(exp(rate))),
           tpm = exp(rate - denom + log(1e6)),
           tpm = replace_na(tpm, 0)) %>%
    ungroup() %>%
    select(sampleid, gene_id, gene_name, tx_id, true_counts, true_tpm = tpm)

## gene-level true expression
ground_truth_genes <- ground_truth_tx_tpms %>%
    group_by(sampleid, gene_id, gene_name) %>%
    summarise(true_counts = sum(true_counts),
              true_tpm = sum(true_tpm)) %>%
    ungroup()

# Compute weighted distances
dist_df <- ground_truth_tx_tpms %>%
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
    summarise(counts = sum(NumReads),
              tpm = sum(TPM)) %>%
    ungroup() %>%
    full_join(ground_truth_genes, by = c("sampleid", "gene_id", "gene_name")) %>%
    mutate_at(vars(counts, true_counts, tpm, true_tpm), ~replace_na(., 0))

## Salmon hla-mapper
salmon_mapper <- sprintf("./pipeline_results/salmon/%s/quant.sf", samples) %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "sampleid") %>%
    left_join(annots, by = c("Name" = "tx_id")) %>%
    group_by(sampleid, gene_id, gene_name) %>%
    summarise(counts = sum(NumReads),
              tpm = sum(TPM)) %>%
    ungroup() %>%
    full_join(ground_truth_genes, by = c("sampleid", "gene_id", "gene_name")) %>%
    mutate_at(vars(counts, true_counts, tpm, true_tpm), ~replace_na(., 0))

## Salmon pers
salmon_pers <- sprintf("./salmon-pers/quant/%s/quant.sf", samples) %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "sampleid") %>%
    left_join(all_annots, by = c("Name" = "tx_id")) %>%
    group_by(sampleid, gene_id, gene_name) %>%
    summarise(counts = sum(NumReads),
              tpm = sum(TPM)) %>%
    ungroup() %>%
    full_join(ground_truth_genes, by = c("sampleid", "gene_id", "gene_name")) %>%
    mutate_at(vars(counts, true_counts, tpm, true_tpm), ~replace_na(., 0))

quant_df <- 
    bind_rows("Salmon ref genome" = salmon_ref,
	      "Salmon personalized" = salmon_pers,
	      "Salmon hla-mapper" = salmon_mapper,
	      .id = "method") %>%
    mutate(method = factor(method, levels = unique(method)))

hla_counts <- quant_df %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    select(method, sampleid, gene_name, counts) %>%
    arrange(method, sampleid, gene_name)

write_rds(hla_counts, "./plot_data/hla_est_counts.rds")

# Assess accuracy

## Compute TPM/TRUE for HLA
count_rates <- quant_df %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    mutate(rate_counts = counts/true_counts,
           rate_tpm = tpm/true_tpm) %>%
    select(method, sampleid, gene_id, gene_name, rate_counts, rate_tpm) %>%
    arrange(method, sampleid, gene_name) %>%
    mutate(gene_name = factor(gene_name, levels = unique(gene_name)))

write_rds(count_rates, "./plot_data/count_rates.rds")








# Fold change in Ref approach
hla_bed <- read_tsv("../analysis/plot_data/hla.bed", col_names = FALSE)

mhc_genes <- read_tsv("../indices/transcript_annotation_df.tsv") %>%
    filter(chr == "chr6") %>%
    mutate(tss = ifelse(strand == "+", start, end)) %>%
    filter(between(tss, min(hla_bed$X2) - 5e5, max(hla_bed$X3) + 5e5)) %>%
    distinct(gene_id, gene_name)

salmon_ref %>%
    inner_join(mhc_genes) %>%
    filter(!(counts < 5 & true_counts < 5)) %>%
    filter(counts > 10) %>%
    mutate(d = counts - true_counts) %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    arrange(sampleid)
    

### HLApers
#hlapers_annot <- sprintf("./hlapers/genotypes/%s_genotypes.tsv", samples) %>%
#    setNames(samples) %>%
#    map_df(read_tsv, .id = "sampleid") %>%
#    mutate(allele = sub("^([^-]+).*$", "\\1", allele),
#	   gene_name = paste0("HLA-", locus)) %>%
#    left_join(distinct(annots, gene_id, gene_name), by = "gene_name") %>%
#    distinct(gene_id, gene_name, tx_id = allele)
#
#hlapers <- sprintf("./hlapers/quant/%s_quant/quant.sf", samples) %>%
#    setNames(samples) %>% .[1] %>%
#    map_df(read_tsv, .id = "sampleid") %>%
#    mutate(Name = sub("^([^-]+).*$", "\\1", Name)) %>%
#    left_join(bind_rows(annots, hlapers_annot), by = c("Name" = "tx_id")) %>%
#    group_by(sampleid, gene_id, gene_name) %>%
#    summarise(counts = sum(NumReads),
#              tpm = sum(TPM)) %>%
#    ungroup() %>%
#    full_join(ground_truth_genes, by = c("sampleid", "gene_id", "gene_name")) %>%
#    mutate_at(vars(counts, true_counts, tpm, true_tpm), ~replace_na(., 0))
#
#
#all_ref_isoforms <- readDNAStringSet("../indices/gencode.transcripts.fa")
#
#hla_all_iso <- read_tsv("../indices/transcript_annotation_df.tsv") %>%
#    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
#    select(gene_name, tx_id)
#
#
#    all_ref_isoforms[hla_all_iso$tx_id] %>%
#    {tibble(tx_id = names(.), len = width(.))} %>%
#    left_join(hla_all_iso) %>%
#    group_by(gene_name) %>%
#    summarise(median_len = median(len)) %>%
#    ungroup()
#      
#
#
#
#hlapers %>%
#    mutate(counts_back = TPM/1e6 * EffectiveLength * sum(NumReads/EffectiveLength)) %>%
#    filter(grepl("HLA-(A|B|C)", gene_name))
#
#
#
#
### hla-mapper::rna
#hlamapper <- sprintf("./hla-mapper/quant/%s_gene.txt", samples) %>%
#    setNames(samples) %>%
#    map_df(~read_tsv(., comment = "#", col_types = "c----dd") %>% 
#	   select(gene_id = 1, len = 2, counts = 3), 
#           .id = "sampleid") %>%
#    left_join(distinct(annots, gene_id, gene_name), by = "gene_id") %>%
#    full_join(ground_truth_genes, by = c("sampleid", "gene_id", "gene_name")) %>%
#    mutate_at(vars(counts:true_counts), ~replace_na(., 0)) %>%
#    select(sampleid, gene_id, gene_name, counts, true_counts)
#
#hlamapper_edit <- sprintf("./hla-mapper/quant/%s_gene_edit.txt", samples) %>%
#    setNames(samples) %>%
#    map_df(~read_tsv(., comment = "#", col_types = "c----dd") %>% 
#	   select(gene_id = 1, len = 2, counts = 3), 
#           .id = "sampleid") %>%
#    left_join(distinct(annots, gene_id, gene_name), by = "gene_id") %>%
#    full_join(ground_truth_genes, by = c("sampleid", "gene_id", "gene_name")) %>%
#    mutate_at(vars(counts:true_counts), ~replace_na(., 0)) %>%
#    select(sampleid, gene_id, gene_name, counts, true_counts)
#
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
