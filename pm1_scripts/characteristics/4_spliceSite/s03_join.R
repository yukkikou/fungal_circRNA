#!/usr/bin/Rscipt

library(tidyverse)

sm = c("C1","C2","EC1","EC3","Y2","Y3","EY3","EY4", "All")
sum_stat = data.frame(matrix(ncol=3)) %>% `colnames<-`(c("Count", "Motif", "Sample"))

for (s in sm){
	tmp = read.table(paste0("stat/",s,".stat.txt")) %>%
		mutate("Sample" = s) %>%
		`colnames<-`(c("Count", "Motif", "Sample"))
	assign(paste0(s,"_stat"), tmp)
	print(tmp)
	sum_stat = rbind(sum_stat, tmp)
}


write_tsv(sum_stat[-1,],"stat/sum_stat.tsv")
save(sum_stat, file = "stat/sum_stat.RData")
