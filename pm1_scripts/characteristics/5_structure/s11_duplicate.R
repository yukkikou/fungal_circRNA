#!/usr/bin/Rscript

library(tidyverse)

circ = read_tsv("circ_structure_uniq.hyphen.tsv", col_names = F)
circ = circ[!duplicated(circ$X1),]

write_tsv(circ, "circ_structure_uniq.hyphen.onetrans.tsv", col_names = F)
