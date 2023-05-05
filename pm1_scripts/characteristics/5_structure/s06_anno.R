#!/usr/bin/Rscript

library(tidyverse)

circ_gene_map = read_tsv("/media/share/node10/disk3/ylab/hxy/PM1_circ_nanopore_202101/5_BSJsite/1_quantify/3_All_l2/stat/2_annotation/circ_gene.sorted.hyphen.map", col_names = c("circ_id","gene_id", "circtrans_id"))

circ_bed13 = read_tsv("circlong_suppIN.gtf.bed13", col_names = c("chr","ss","es","circ_id","strand",
                                                    "block_number","block_size","block_ss","block_es"))

circ_merge = circ_bed13 %>%
        left_join(circ_gene_map, by = c("circ_id" = "circ_id"))

circ_attr = data.frame(
    circ_attr_circ_id = paste0("circ_id ","|", circ_merge$circ_id, '|;'),
    circ_attr_gene_id = paste0("gene_id ", "|",circ_merge$gene_id, '|;'),
    circ_attr_circtrans_id = paste0("circtrans_id ", "|",circ_merge$circtrans_id, '|;')
)

for (i in 1:dim(circ_attr)[1]){
    circ_attr[i, "str"] = paste(circ_attr[i,"circ_attr_circ_id"],circ_attr[i,"circ_attr_gene_id"],circ_attr[i,"circ_attr_circtrans_id"])
}

print(head(circ_attr))

circ_gtf = data.frame(
    chr = circ_merge$chr,
    src = ".",
    feat = "circExon",
    ss = as.numeric(circ_merge$ss) + 1,
    es = as.numeric(circ_merge$es),
    score = ".",
    strand = circ_merge$strand,
    phase = ".",
    attr = circ_attr$str
)

write_tsv(circ_gtf, "circlong_suppIN.final.sorted.gtf", col_names = F)
