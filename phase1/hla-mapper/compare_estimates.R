library(tidyverse)

star <- read_tsv("./quant/STAR/66K00003/quant.sf")
xlate <- read_tsv("./quant/xlate/66K00003/quant.sf")

expression_df <- 
    left_join(select(star, Name, star = TPM),
              select(xlate, Name, xlate = TPM))
    
expression_df %>%
    mutate_at(vars(star, xlate), ~`+`(., 1/1e5)) %>%
    ggplot(aes(xlate, star)) +
    geom_jitter(alpha = .1) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() +
    labs(x = expression(paste("log"[10], italic(TPM ("sam-xlate")))),
         y = expression(paste("log"[10], italic(TPM ("STAR-toTranscriptome")))))


ggsave("estimates.png", height = 5, width = 5)

expression_df %>%
    filter(star > 0, xlate == 0)

expression_df %>%
    filter(star > 1, xlate == 0)

expression_df %>%
    filter(xlate > 0, star == 0)

expression_df %>%
    filter(xlate > 1, star == 0)


gtf <- "/home/vitor/gencode_data/v36/gencode.v36.annotation.gtf.gz" %>%
    read_tsv(comment = "##", col_names = FALSE) %>%
    setNames(c("seqname", "source", "feature", "start", "end", 
               "score", "strand", "frame", "attribute"))

gene_transc_map <- gtf %>%
    filter(feature == "transcript") %>%
    transmute(gene_id = str_extract(attribute, "(?<=gene_id\\s\")[^\"]+"),
              transcript_id = str_extract(attribute, "(?<=transcript_id\\s\")[^\"]+"))

expression_genes <- expression_df %>%
    left_join(gene_transc_map, by = c("Name" = "transcript_id")) %>%
    group_by(gene_id) %>%
    summarise(star = sum(star),
              xlate = sum(xlate)) %>%
    ungroup()

expression_genes %>%
    mutate_at(vars(star, xlate), ~`+`(., 1/1e5)) %>%
    ggplot(aes(xlate, star)) +
    geom_jitter(alpha = .1) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() +
    labs(x = expression(paste("log"[10], italic(TPM ("sam-xlate")))),
         y = expression(paste("log"[10], italic(TPM ("STAR-toTranscriptome")))))

