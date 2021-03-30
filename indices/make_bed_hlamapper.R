library(tidyverse)

bedstar <- "./STAR/chrNameLength.txt" %>%
    read_tsv(col_names = FALSE) %>%
    filter(X1 != 6 & X1 != "chr6") %>%
    mutate(start = 1L) %>%
    select(chr = X1, start, end = X2)

bedhlam <- 
    "~/HLAMAPPER/hlamapper/hla-mapper_db_004.1_HLA/bed/not_target_rna.bed" %>%
    read_tsv(col_names = FALSE) %>%
    filter(X1 == 6 | X1 == "chr6") %>%
    select(chr = X1, start = X2, end = X3)

final_bed <- bind_rows(bedhlam, bedstar)

write_tsv(final_bed,
	  "~/HLAMAPPER/hlamapper/hla-mapper_db_004.1_HLA/bed/not_target_rna.bed",
	  col_names = FALSE)
