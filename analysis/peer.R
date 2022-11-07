library(tidyverse)
library(peer)

phenotypes <- 
    read_tsv("./salmon-pers/quants.bed.gz") %>%
    select(gid, starts_with("66K")) %>%
    pivot_longer(-gid, names_to = "sampleid", values_to = "tpm") %>%
    group_by(gid) %>%
    filter(mean(tpm > 1) > .75) %>%
    ungroup() %>%
    pivot_wider(names_from = gid, values_from = "tpm") %>%
    column_to_rownames("sampleid") %>%
    data.matrix()

covar_matrix <- 
    tibble(sampleid = rownames(phenotypes)) %>%
    mutate(timepoint = case_when(grepl("t1$", sampleid) ~ 0,
				 grepl("t2$", sampleid) ~ 1)) %>%
    as.data.frame() %>%
    column_to_rownames("sampleid") %>%
    data.matrix()

model <- PEER()
PEER_setPhenoMean(model, phenotypes)
PEER_setCovariates(model, covar_matrix)
PEER_setNk(model, 25)
PEER_setNmax_iterations(model, 1e5)
PEER_setAdd_mean(model, TRUE)
PEER_update(model)

#cor(PEER_getX(model)[,1], PEER_getCovariates(model)[,1])

png("./peer_model.png")
PEER_plotModel(model)
dev.off()

resids <- PEER_getResiduals(model)
dimnames(resids) <- dimnames(phenotypes)

resid_df <- resids %>%
    as.data.frame() %>%
    rownames_to_column("sampleid") %>%
    as_tibble() %>%
    pivot_longer(-sampleid, names_to = "gene_id", values_to = "resid")

write_tsv(resid_df, "./plot_data/peer_residuals.tsv")

factors <- PEER_getX(model)
rownames(factors) <- rownames(phenotypes)
as_tibble(factors, rownames = "sampleid") %>%
    write_tsv("./plot_data/peer_factors.tsv")
