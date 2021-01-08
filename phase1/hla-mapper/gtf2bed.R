library(tidyverse)

gtf <- "/home/vitor/gencode_data/v36/gencode.v36.annotation.gtf.gz" %>%
    read_tsv(comment = "##", col_names = FALSE) %>%
    setNames(c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

tmp <- gtf %>%
    filter(feature == "exon") %>%
    mutate(start = start - 1L,
	   transcript_id = str_extract(attribute, "(?<=transcript_id\\s\")[^\"]+"),
	   block_size = end - start) %>%
    select(seqname, transcript_id, start, end, strand, block_size) %>%
    mutate(seqname = factor(seqname, levels = paste0("chr", c(1:22, "X", "Y", "M")))) %>%
    arrange(seqname, start) %>%

bed <- tmp %>%
    group_by(seqname, transcript_id) %>%
    mutate(block_start = start - first(start)) %>%
    summarise(start = first(start),
	      end = last(end),
	      strand = first(strand),
	      n = n(),
	      block_sizes = paste(block_size, collapse = ","),
	      block_starts = paste(block_start, collapse = ",")) %>%
    ungroup() %>%
    mutate(across(block_sizes:block_starts, ~paste0(., ","))) %>%
    arrange(seqname, start) %>%
    mutate(score = 0, thickStart = 0, thickEnd = 0, itemRgb = 0) %>%
    select(chrom = seqname, start, end, name = transcript_id, 
	   score, strand, thickStart, thickEnd, itemRgb, 
	   blockCount = n, blockSizes = block_sizes, blockStarts = block_starts)

write_tsv(bed, "./gencode_v36.bed", col_names = FALSE)
