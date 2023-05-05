#!/usr/bin/Rscript

library(tidyverse)

circ_bed = read_tsv("circlong_suppIN.r.bed12", col_names = c("chr","ss","es","circ_id","score","strand","x7","x8","x9","block_num","block_size","block_ss"))

single_id = read_tsv("unAnnoSingle.id", col_names = "circ_id")
aban_id = read_tsv("unAnnotation.id", col_names = "circ_id")

circ_bed_clean = circ_bed %>% filter(!(circ_id %in% aban_id$circ_id)) %>% as.data.frame()

for (i in 1:dim(circ_bed_clean)[1]){
    circ_tmp = circ_bed_clean[i,'circ_id']
    circ_bed_clean[i,'x7'] = circ_bed_clean[i,'ss']
    circ_bed_clean[i,'x8'] = circ_bed_clean[i,'es']
    if (circ_tmp %in% single_id$circ_id){
        print(circ_tmp)
        circ_bed_clean[i,'block_ss'] = circ_bed_clean[i,'ss']
        circ_bed_clean[i,'block_size'] = circ_bed_clean[i,'es'] - circ_bed_clean[i,'ss']
        print(circ_bed_clean[i,])
    }
}

write_tsv(circ_bed_clean,"circlong_suppIN.r.clean.bed12", col_names = F)
