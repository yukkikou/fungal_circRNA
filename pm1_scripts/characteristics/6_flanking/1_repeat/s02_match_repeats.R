#!/usr/bin/Rscript

library(tidyverse)

upinter = read_tsv("flk.up.rep.inter", col_names = F)
downinter = read_tsv("flk.down.rep.inter", col_names = F)
cols = c("chr","ss","es","circID","geneID","strand","rchr","rss","res","rID","rtype","rstand","inter")

colnames(upinter) = cols
colnames(downinter) = cols

circ_simple_inter = inner_join(upinter, downinter,
	by = c("circID" = "circID"),
	suffix = c(".up", ".down"))

circ_full_inter = inner_join(upinter, downinter,
	by = c("circID" = "circID", "rID" = "rID"), 
	suffix = c(".up", ".down"))


write_tsv(circ_simple_inter, "circ_simple_inter.tsv")
write_tsv(circ_full_inter, "circ_full_inter.tsv")
