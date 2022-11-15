library(Biostrings)
library(tidyverse)

gtffile <- "/raid/genevol/gencode/gencode.v37.primary_assembly.annotation.gtf"

annots <- read_tsv(gtffile, comment = "#", col_names = FALSE, col_types = "ccciicccc")

gene_annots <- annots %>%
    filter(X3 == "gene") %>%
    filter(grepl("HLA-A", X9))

hla_gtf <- annots %>%
    filter(grepl("gene_name \"HLA-A\";", X9))

# write GTF header
read_lines(gtffile, n_max = 5) %>%
    write_lines("./hla.gtf")

# write GTF contents
write_tsv(hla_gtf, "./hla.gtf", append = TRUE, escape = "none")

# mask genome
genome38 <- readDNAStringSet("/raid/genevol/gencode/GRCh38.primary_assembly.genome.fa")

chr6 <- genome38["chr6 6"] 

proximal <- subseq(chr6, 1, min(gene_annots$X4) - 1000)
terminal <- subseq(chr6, max(gene_annots$X5) + 1000, width(chr6))

subseq(chr6, 1, min(gene_annots$X4) - 1000) <- 
   DNAString(paste(rep("N", width(proximal)), collapse = ""))

subseq(chr6, max(gene_annots$X5) + 1000, width(chr6)) <- 
    DNAString(paste(rep("N", width(terminal)), collapse = ""))

writeXStringSet(chr6, "chr6_masked.fasta")
