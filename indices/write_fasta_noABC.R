library(Biostrings)
library(tidyverse)

annots <- read_tsv("../../indices/transcript_annotation_df.tsv") %>%
    filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C"))

index <- readDNAStringSet("../../indices/gencode.transcripts.fa") %>%
    .[! names(.) %in% annots$tx_id ]

writeXStringSet(index, "transcripts_noABC.fa")
