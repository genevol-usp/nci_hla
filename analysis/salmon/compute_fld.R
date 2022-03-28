library(tidyverse)

process_fld <- function(x) {
    fld <- scan(x)
    
    tibble(mu = sum(1:length(fld) *  fld)) %>%
	mutate(sigma = sqrt(sum( ( (1:length(fld) - mu)^2 ) * fld )))
}

ids <- dir("quant")

fld <- sprintf("./quant/%s/libParams/flenDist.txt", ids) %>%
    setNames(ids) %>%
    map_df(process_fld, .id = "sampleid") %>%
    separate(sampleid, c("sampleid", "timepoint"), sep = "_")


fld %>%
    group_by(timepoint) %>%
    summarise_at(vars(mu, sigma), mean)
