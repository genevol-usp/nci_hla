library(Biostrings)
library(tidyverse)

# sample IDs
samples <- readLines("./samples.txt")

# compute sequence divergence of personalized transcripts
# in regard to original sequences
transcripts <- "/home/vitor/simulation_rnaseq/simulation_index.fa" %>%
    readDNAStringSet()

pers_isoforms <- transcripts[grepl("\\*", names(transcripts))] 

hla_transcripts <- names(pers_isoforms) %>%
    str_split("_") %>%
    map(1) %>%
    unlist() %>%
    unique()

ref_isoforms <- "./salmon-reads/ensembl.transcripts.fa" %>%
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
annots <- "/home/vitor/hisat2/grch38_snp_tran/Homo_sapiens.GRCh38.99.gtf" %>% 
    read_tsv(comment = "#", col_names = FALSE, col_types = "c-cii-c-c") %>%
    filter(X3 == "transcript") %>%
    transmute(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
              gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
              tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"))

# True counts

lengths_df <- tibble(tx_id = names(transcripts),
		     len = width(transcripts),
		     eff_len = len - 300L + 1L)

ground_truth <- read_tsv("/home/vitor/simulation_rnaseq/phenotypes.tsv") %>%
    add_column(tx_id = names(transcripts)) %>%
    left_join(lengths_df, by = "tx_id") %>%
    pivot_longer(-(tx_id:eff_len), names_to = "sampleid", values_to = "true_counts") %>% 
    mutate(rate = true_counts/eff_len,
	   rate = replace_na(rate, 0),
	   rate = if_else(eff_len == 0 | eff_len < 0, 0, rate)) %>% 
    group_by(sampleid) %>%
    mutate(tpm = rate/sum(rate) * 1e6) %>%
    ungroup()

true_hla_tpm <- ground_truth %>% 
    filter(grepl("\\*", tx_id)) %>%
    separate(tx_id, c("tx_id", "allele"), sep = "_") %>%
    mutate(gene_name = sub("^([^*]+)\\*.+$", "HLA-\\1", allele)) %>%
    group_by(sampleid, gene_name) %>%
    summarise(true_tpm = sum(tpm)) %>%
    ungroup()


# Multiple Salmon quantification methods
salmon_nocorr <- file.path("./salmon-reads/quant_nocorrection", samples, "quant.sf") %>%
    setNames(sprintf("salmon-reads_NoBiasCorrect.%s", samples))

salmon_corr <- file.path("./salmon-reads/quant_biascorrection", samples, "quant.sf") %>%
    setNames(sprintf("salmon-reads_BiasCorrect.%s", samples))

salmon_aln_nocorr <- file.path("./salmon-alignment/quant_nocorrection", samples, "quant.sf") %>%
    setNames(sprintf("salmon-aln_NoBiasCorrect.%s", samples))
 
salmon_aln_corr <- file.path("./salmon-alignment/quant_biascorrection", samples, "quant.sf") %>%
    setNames(sprintf("salmon-aln_BiasCorrect.%s", samples))

salmon_df <- c(salmon_nocorr, salmon_corr, salmon_aln_nocorr, salmon_aln_corr) %>%
    map_df(read_tsv, .id = "ID") %>%
    extract(ID, c("method", "sampleid"), "(.+)\\.(.+)")

# compute MARDs

## First, standardize HLA transcripts
ground_truth_std <- ground_truth %>%
    mutate(tx_id = sub("(_.+)$", "", tx_id)) %>%
    group_by(sampleid, tx_id) %>%
    summarise(true = sum(tpm)) %>%
    ungroup()

## compute MARD
ard_df <- salmon_df %>%
    select(method, sampleid, tx_id = Name, tpm = TPM) %>%
    full_join(ground_truth_std, by = c("sampleid", "tx_id")) %>%
    mutate_at(vars(tpm, true), ~replace_na(., 0)) %>%
    mutate(ard = if_else(tpm == 0 & true == 0, 0, abs(true - tpm) / (true + tpm)))

mard_genomewide <- ard_df %>%
    group_by(method, sampleid) %>%
    summarise(mard = mean(ard)) %>%
    ungroup()

mard_hla <- ard_df %>%
    filter(tx_id %in% hla_transcripts$tx_id) %>%
    group_by(method, sampleid) %>%
    summarise(mard = mean(ard)) %>%
    ungroup()

write_rds(mard_genomewide, "./plot_data/mard_genomewide.rds")
write_rds(mard_hla, "./plot_data/mard_hla.rds")




