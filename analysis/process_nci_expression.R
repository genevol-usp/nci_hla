library(tidyverse)
library(readxl)

nci_genos <- "/raid/genevol/nci_rnaseq/genotypes_nci_tidy.tsv" %>%
    read_tsv() %>%
    mutate(allele = sub("^([ABC]\\*)(\\d+:)(\\d+).*$", "\\1\\2\\3", allele),
           lineage = sub("^([^:]+).*$", "\\1", allele),
           gene_name = sub("HLA-", "", gene_name)) %>%
    unite(allele_i, gene_name:i, sep = "") %>%
    select(-allele) %>%
    pivot_wider(names_from = allele_i, values_from = lineage) %>%
    arrange(sampleid)

nci_expression <- 
    "/raid/genevol/nci_rnaseq/phase1/HLA-A, -B, -C expression levels for Pat.xlsx" %>%
    read_excel() %>%
    janitor::clean_names() %>%
    select(sampleid = 3, ends_with("rna"))

nci_df <- inner_join(nci_genos, nci_expression, by = "sampleid")

# HLA-A
a_expr <- select(nci_df, sampleid, A1, A2, rna = hla_a_m_rna) %>%
    drop_na()

## double file for allele symmetry
a_expr_1 <- a_expr
a_expr_2 <- a_expr

a_expr_2$A1 <- a_expr$A2
a_expr_2$A2 <- a_expr$A1

a_expr_dbl <- rbind(a_expr_1, a_expr_2) %>%
    arrange(sampleid)

## Linear models
allele_effect <- lm(rna ~ A1 + A2, data = a_expr_dbl)
summary(allele_effect)

## Get predictions for models
pred_level <- predict(allele_effect, se.fit = T)
a_expr_dbl$predicted <- pred_level$fit
a_expr_dbl$pred_se <- pred_level$se.fit

x_df <- a_expr_dbl %>% 
    filter(predicted < 0) %>%
    arrange(predicted)

## adhoc correction
x <- (x_df$rna[1] - x_df$predicted[1]) / x_df$pred_se[1]
a_expr_dbl$predicted_adj <- pred_level$fit + pred_level$se.fit*x




# HLA-B
b_expr <- select(nci_df, sampleid, B1, B2, rna = hla_b_m_rna) %>%
    drop_na()

## double file for allele symmetry
b_expr_1 <- b_expr
b_expr_2 <- b_expr

b_expr_2$B1 <- b_expr$B2
b_expr_2$B2 <- b_expr$B1

b_expr_dbl <- rbind(b_expr_1, b_expr_2) %>%
    arrange(sampleid)

## Linear models
allele_effect_b <- lm(rna ~ B1 + B2, data = b_expr_dbl)
summary(allele_effect_b)

## Get predictions for models
pred_level_b <- predict(allele_effect_b, se.fit = T)
b_expr_dbl$predicted <- pred_level_b$fit
b_expr_dbl$pred_se <- pred_level_b$se.fit

b_expr_dbl %>% 
    filter(predicted < 0)

## No adhoc correction needed




# HLA-C
c_expr <- select(nci_df, sampleid, C1, C2, rna = hla_c_m_rna) %>%
    mutate(rna = as.numeric(rna)) %>%
    drop_na()

## double file for allele symmetry
c_expr_1 <- c_expr
c_expr_2 <- c_expr

c_expr_2$C1 <- c_expr$C2
c_expr_2$C2 <- c_expr$C1

c_expr_dbl <- rbind(c_expr_1, c_expr_2) %>%
    arrange(sampleid)

## Linear models
allele_effect_c <- lm(rna ~ C1 + C2, data = c_expr_dbl)
summary(allele_effect_c)

## Get predictions for models
pred_level_c <- predict(allele_effect_c, se.fit = T)
c_expr_dbl$predicted <- pred_level_c$fit
c_expr_dbl$pred_se <- pred_level_c$se.fit

c_expr_dbl %>% 
    filter(predicted < 0)

## No adhoc correction needed


# Final dataset

a_final <- a_expr_dbl %>%
    distinct(sampleid, .keep_all = TRUE) %>%
    select(sampleid, A1, A2, rna = predicted_adj) %>%
    pivot_longer(A1:A2, names_to = "gene_name", values_to = "allele") %>%
    mutate(gene_name = paste0("HLA-", gene_name),
           gene_name = sub("\\d$", "", gene_name)) %>%
    select(sampleid, gene_name, allele, rna)

b_final <- b_expr_dbl %>%
    distinct(sampleid, .keep_all = TRUE) %>%
    select(sampleid, B1, B2, rna = predicted) %>%
    pivot_longer(B1:B2, names_to = "gene_name", values_to = "allele") %>%
    mutate(gene_name = paste0("HLA-", gene_name),
           gene_name = sub("\\d$", "", gene_name)) %>%
    select(sampleid, gene_name, allele, rna)
    
c_final <- c_expr_dbl %>%
    distinct(sampleid, .keep_all = TRUE) %>%
    select(sampleid, C1, C2, rna = predicted) %>%
    pivot_longer(C1:C2, names_to = "gene_name", values_to = "allele") %>%
    mutate(gene_name = paste0("HLA-", gene_name),
           gene_name = sub("\\d$", "", gene_name)) %>%
    select(sampleid, gene_name, allele, rna)

final_df <- bind_rows(a_final, b_final, c_final) %>%
    arrange(sampleid, gene_name)

write_rds(final_df, "./plot_data/nci_allele.rds")
