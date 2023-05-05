#!/usr/bin/Rscript

#library(tidyverse)

sample=c(paste0("M",1:4),paste0("Y",1:4))
files=paste("IN", sample, "Quant.ratio", sep = "_")

all_ratio = c()
for (i in files){
	print(i)
	tmp = as.numeric(readLines(i))	
	print(quantile(tmp, probs = seq(0,1,0.1)))
	all_ratio = c(all_ratio, tmp)
}

message("All ratio")
quantile(all_ratio, probs = seq(0,1,0.03))
