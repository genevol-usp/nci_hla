library(tidyverse)

sam1 <- read_tsv("./66K00003_Aligned.sam", col_names = FALSE)
sam2 <- read_tsv("./66K00003.sam", col_names = FALSE)
sam3 <- read_tsv("./66K00003_originalmapper.sam", col_names = FALSE)


sam3$X1[! sam3$X1 %in% sam1$X1]


sam1 %>% count(X2, X16, sort = TRUE) %>% print(n = Inf)

grp <- 
    tribble(~grp, ~flag, ~mate,
	    "1", 83, 1,
	    "1", 163, 2,
	    "2", 338, 1,
	    "2", 419, 2,
	    "3", 99, 1,
	    "3", 147, 2,
	    "4", 355, 1,
	    "4", 403, 2)

sam1 %>% 
    select(read = X1, flag = X2, start = X4, qual = X5, cigar = X6, mate_start = X8, len = X9, X16) %>%
    filter(qual == 255) %>%
    filter(start > 29941000) %>%
    arrange(read)
