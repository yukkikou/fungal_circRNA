#!/usr/bin/Rscript

# option
options(stringsAsFactors=F)

# packages
library(tidyverse)

# function
basePre = function(bases){
	base_table = table(unlist(bases)) |>
		as.data.frame() |>
		`colnames<-`(c("Base", "Count")) |>
		mutate("Precent" = Count/length(bases))
	base_table
}

# body
motif_table = read_tsv("all.motifRev.tsv")
motif_seq = str_split(unlist(motif_table[,10]), pattern="", simplify = T)

apply(motif_seq, 2, basePre)
