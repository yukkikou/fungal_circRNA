#!/usr/bin/Rscript

library(tidyverse)

trans_cpmMean = read_tsv("../../7_expression/2_linear/1_featureCounts/results/transcript_cpmMean.clean.tsv")

trans_cpmMean_out = trans_cpmMean %>% 
    mutate(gene_id = paste0(str_split(trans_cpmMean$transcript_id, pattern = "\\.", simplify = T)[,1], ".00"),
            cpmMax = as.numeric(apply(trans_cpmMean, 1, \(x){max(x[2:3])}, simplify = T)))

trans_keep = trans_cpmMean_out %>% 
    group_by(gene_id) %>%
    summarise(max = max(cpmMax))

trans_final = trans_cpmMean_out %>% inner_join(trans_keep, by = c("gene_id" = "gene_id", "cpmMax" = "max")) %>% 
    filter(cpmMax != 0)

write_tsv(trans_final, "trans_highest.tsv")`
