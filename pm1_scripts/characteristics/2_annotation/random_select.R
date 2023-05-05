#!/usr/bin/Rscript

library(tidyverse)

circ_length = read_tsv("circ_IN.length", col_names = c("circ_id","exon","length"))

idx = sample(1:dim(circ_length)[1], 30)

write_tsv(circ_length[idx,], "random_circ.length")
