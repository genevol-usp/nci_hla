library(tidyverse)

samples <- read_tsv("./samples.txt", col_names = c("timepoint", "sampleid"))
sample_ids <- sprintf("%s_t%s", samples$sampleid, samples$timepoint)

writeLines(sample_ids, "sample_ids.txt")

gene_annots <- read_tsv("../indices/gene_annotation_df.tsv")
tx_annots <- read_tsv("../indices/transcript_annotation_df.tsv")

###############################################################################
# Salmon ref
salmon_ref <- file.path("./salmon/quant", sample_ids[1:96], "quant.sf") %>%
    setNames(sub("_t\\d$", "", sample_ids[1:96])) %>%
    map_df(read_tsv, .id = "sampleid") 

salmon_ref_genes <- salmon_ref %>%
    left_join(tx_annots, by = c("Name" = "tx_id")) %>%
    group_by(sampleid, gene_id, gene_name) %>%
    summarise(tpm = sum(TPM)) %>%
    ungroup()
	   
salmon_ref_bed <- salmon_ref_genes %>%
    pivot_wider(names_from = sampleid, values_from = tpm) %>%
    left_join(gene_annots, by = c("gene_id", "gene_name")) %>%
    mutate(chr = factor(chr, levels = unique(gene_annots$chr))) %>%
    select(`#chr` = chr, start, end, id = gene_name, gid = gene_id, strd = strand, starts_with("66K")) %>%
    arrange(`#chr`, start)


write_tsv(salmon_ref_bed, "./salmon/quants.bed")

sub("_t\\d$", "", sample_ids[1:96]) %>%
    write_lines("./sample_ids_t1.txt")


salmon_ref_tx <- salmon_ref %>%
    group_by(tx_id = Name) %>%
    summarise(tpm = mean(TPM),
	      efflen = mean(EffectiveLength)) %>%
    ungroup()

write_rds(salmon_ref_tx, "./plot_data/salmon_ref_tx_summary.rds")



###############################################################################
# Salmon Pers 
salmon_pers <- file.path("./salmon-pers/quant", sample_ids, "quant.sf") %>%
    setNames(sample_ids) %>%
    map_df(read_tsv, .id = "sampleid") %>%
    mutate(tx_id = ifelse(grepl("\\_[ABC]\\*", Name), 
                          sub("^([^_]+).*$", "\\1", Name),
                          Name)) %>%
    left_join(tx_annots, by = "tx_id") %>%
    select(sampleid, Name, tx_id, gene_id, gene_name, EffectiveLength, NumReads, TPM)

salmon_pers_gene <- salmon_pers %>% 
    group_by(sampleid, gene_id, gene_name) %>%
    summarise(tpm = sum(TPM)) %>%
    ungroup()

salmon_pers_bed <- salmon_pers_gene %>%
    pivot_wider(names_from = sampleid, values_from = tpm) %>%
    left_join(gene_annots, by = c("gene_id", "gene_name")) %>%
    mutate(chr = factor(chr, levels = unique(gene_annots$chr))) %>%
    select(`#chr` = chr, start, end, id = gene_name, gid = gene_id, strd = strand, 
	   starts_with("66K")) %>%
    arrange(`#chr`, start)

write_tsv(salmon_pers_bed, "./salmon-pers/quants.bed")

covars <- tibble(ids = sample_ids) %>%
    mutate(batch = sub("^66K\\d+_(t[12])$", "\\1", ids)) %>%
    arrange(ids) %>%
    pivot_longer(-ids, names_to = "id") %>%
    pivot_wider(names_from = ids, values_from = value)

write_tsv(covars, "./salmon-pers/technical_covars.txt")

salmon_pers %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    write_tsv("./plot_data/salmon_pers_hla_isoforms.tsv")


# B2M normalization
salmon_pers_b2m <- salmon_pers %>% 
    filter(grepl("_t1", sampleid)) %>%
    mutate(sampleid = sub("_t\\d$", "", sampleid),
	   tx_id = ifelse(grepl("\\_[ABC]\\*", Name), 
			  sub("^([^_]+).*$", "\\1", Name),
			  Name)) %>%
    left_join(tx_annots, by = "tx_id") %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C", "B2M")) %>%
    group_by(sampleid, gene_name) %>%
    summarise(counts = sum(NumReads)) %>%
    ungroup() %>%
    pivot_wider(names_from = gene_name, values_from = counts) %>%
    pivot_longer(contains("HLA"), names_to = "gene_name", values_to = "counts") %>%
    group_by(gene_name) %>%
    mutate(norm_counts = scale(counts/B2M)[,1]) %>%
    ungroup() %>%
    select(sampleid, gene_name, norm_counts)

write_rds(salmon_pers_b2m, "./plot_data/salmon_pers_b2m.rds")

# HLA allele-level expression
salmon_allele <- salmon_pers %>%
    group_by(sampleid) %>%
    mutate(cpm = NumReads/sum(NumReads) * 1e6) %>%
    ungroup() %>%
    filter(grepl("\\_[ABC]\\*", Name)) %>%
    extract(Name, "allele", ".+_(.+)") %>%
    group_by(sampleid, allele) %>%
    summarise(cpm = sum(cpm),
              tpm = sum(TPM)) %>%
    ungroup() %>%
    mutate(gene_name = sub("^([^*]+).+$", "HLA-\\1", allele)) %>%
    group_by(sampleid, gene_name) %>%
    mutate(n = n_distinct(allele),
           zyg = case_when(n == 1L ~ "1_1",
                           n == 2L ~ "1",
                           TRUE ~ NA_character_)) %>%
    ungroup() %>%
    mutate(cpm = ifelse(zyg == "1_1", cpm/2L, cpm),
           tpm = ifelse(zyg == "1_1", tpm/2L, tpm)) %>%
    separate_rows(zyg, sep = "_") %>%
    select(sampleid, gene_name, allele, cpm, tpm)

write_rds(salmon_allele, "./plot_data/salmon_pers_allele.rds")

# HLApers
file.path("./hlapers/quant", sample_ids, "quant.sf") %>%
    setNames(sample_ids) %>% .[1] %>%
    map_df(read_tsv, .id = "sampleid")





###############################################################################
# hla-mapper
hlamapper_all <- "./hla-mapper/quant/%s_gene_rmPseudo.txt" %>%
    sprintf(sample_ids[1:96]) %>%
    setNames(sub("_t\\d$", "", sample_ids[1:96])) %>%
    map_df(~read_tsv(., comment = "#") %>%
	   select(gene_id = Geneid, Length, counts = contains("66K")) %>%
	   mutate(rate = log(counts) - log(Length),
		  denom = log(sum(exp(rate))),
		  tpm = exp(rate - denom + log(1e6))) %>%
	   select(gene_id, Length, counts, tpm),
           .id = "sampleid")
    
hlamapper_bed <- hlamapper_all %>%
    select(sampleid, gene_id, tpm) %>%
    pivot_wider(names_from = sampleid, values_from = tpm) %>%
    left_join(gene_annots, by = "gene_id") %>%
    mutate(chr = factor(chr, levels = unique(gene_annots$chr))) %>%
    select(`#chr` = chr, start, end, id = gene_name, gid = gene_id, strd = strand, starts_with("66K")) %>%
    arrange(`#chr`, start)

write_tsv(hlamapper_bed, "./hla-mapper/quants.bed")

hlamapper <- hlamapper_all %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    select(sampleid, gene_name, counts)

# Exon-based quants
hlamatch <- read_tsv("../indices/personalize_transcripts/genome_match.tsv") %>%
    mutate(pos = case_when(locus == "A" ~ map2(start, end, ~.x:.y),
                           locus == "B" | locus == "C" ~ map2(start, end, ~.y:.x))) %>%
    select(locus, allele, pos) %>%
    unnest(pos) %>%
    group_by(locus, allele) %>%
    mutate(i = seq_len(n())) %>%
    ungroup()

hla_coords <- unique(hlamatch$locus) %>%
    map_df(~hlaseqlib::hla_read_alignment(., imgtdb = "~/IMGTHLA", 
                                          imgtfile = "gen", by_exon = TRUE) %>%
               filter(allele %in% unique(hlamatch$allele)) %>%
               mutate(cds = gsub("\\.|\\*", "", cds),
                      cds = str_split(cds, "")) %>%
               unnest(cds)) %>%
    group_by(allele) %>%
    mutate(i = seq_len(n())) %>%
    ungroup()

hla_exon_coords <- left_join(hla_coords, hlamatch, by = c("allele", "i")) %>%
    mutate(gene_name = paste0("HLA-", locus)) %>%
    unite(feature, c("feature", "idx_grp"), sep = " ") %>%
    select(gene_name, feature, i, pos) %>%
    filter(feature %in% c("exon 2", "exon 3")) %>%
    group_by(gene_name, feature) %>%
    summarise(start = min(pos), end = max(pos)) %>%
    ungroup()
 
exon_quants <- "./hla-mapper/quant/%s_exon.txt" %>%
    sprintf(sample_ids[1:96]) %>%
    setNames(sample_ids[1:96]) %>%
    map_df(. %>% 
               read_tsv(comment = "#") %>%
               rename("counts" = 7) %>%
               filter(Chr == "chr6") %>%
               setNames(tolower(names(.))) %>%
               inner_join(hla_exon_coords, by = c("start", "end")),
           .id = "sampleid")

exon_normd <- exon_quants %>%
    group_by(sampleid, gene_name) %>%
    summarise(counts = weighted.mean(counts, length)) %>%
    ungroup() %>%
    mutate(sampleid = sub("_t\\d$", "", sampleid)) %>%
    left_join(filter(hlamapper, gene_name == "B2M") %>% select(sampleid, b2m = counts)) %>%
    mutate(b2m_norm_counts = counts/b2m) %>%
    select(sampleid, gene_name, counts, b2m_norm_counts)

hlamapper %>%
    rename(gene_counts = counts) %>%
    left_join(exon_normd) %>%
    filter(gene_name != "B2M") %>%
    rename(exon_counts = counts, b2m_norm_exon_counts = b2m_norm_counts) %>%
    write_rds("./plot_data/hlamapper.rds")



