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
    
    hla_read_alignment(hla_locus, IMGTDB, imgtfile = "gen", exons = 2:3, by_exon = TRUE) %>%
	mutate(cds = map_chr(cds, replace_dots),
	       cds = gsub("\\*", "+", cds),
	       cds = gsub("\\.", "-", cds))
}


# HLA loci
loci <- c("A", "B", "C")

IMGTDB <- "/home/vitor/IMGTHLA"

# Gene annotations
annots <- "/raid/genevol/gencode/gencode.v37.primary_assembly.annotation.gtf" %>% 
    read_tsv(comment = "#", col_names = FALSE, col_types = "c-cii-c-c")

gene_annots <- annots %>%
    filter(X3 == "gene") %>%
    transmute(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
              gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
              strand = X7) %>%
    filter(gene_name %in% paste0("HLA-", loci))

# match to genome
ref_alleles <- read_tsv("../indices/personalize_transcripts/genome_match.tsv")

hlagen <- map_df(loci, read_hla) %>% 
    filter(allele %in% ref_alleles$allele) %>%
    mutate(gene_name = sub("^([^*]+).+$", "HLA-\\1", allele)) %>%
    unite(allele, c(allele, idx_grp), sep = "_") %>%
    left_join(select(gene_annots, gene_name, strand), by = "gene_name") %>%
    select(gene_name, strand, allele, cds)

hlagen$cds[hlagen$strand == "-"] <- 
    as.character(reverseComplement(DNAStringSet(hlagen$cds[hlagen$strand == "-"])))

hlaset <- gsub("[+-]", "", hlagen$cds) %>%
    DNAStringSet() %>%
    setNames(hlagen$allele)

genome <- readDNAStringSet("/raid/genevol/gencode/GRCh38.primary_assembly.genome.fa")

mhc <- subseq(genome["chr6 6"], 29e6, 32e6)

genomematch <- matchPDict(hlaset, mhc[[1]]) %>%
    as.list() %>%
    map_df(as.data.frame, .id = "allele") %>%
    as_tibble() %>%
    mutate_at(vars(start, end), ~`+`(., 29e6L - 1L)) %>%
    mutate(locus = sub("^([^*]+).+$", "\\1", allele)) %>%
    select(locus, allele, start, end, width)

genomematch %>%
    mutate(gene_name = paste0("HLA-", locus)) %>%
    separate(allele, c("allele", "exon"), sep = "_") %>%
    select(gene_name, exon, start, end) %>%
    write_tsv("./plot_data/hla_exon_coords.tsv")



