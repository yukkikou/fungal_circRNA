# packages
library(dplyr)
library(purrr)
library(magrittr)
library(stringr)
library(ggplot2)
library(clusterProfiler)
library(RColorBrewer)
library(DOSE)

pm1GeneGo = read_tsv("D://_CircRNA/scripts/Enrichment/pm1.go")
pm1ReportedGene = read_tsv("D://_CircRNA/scripts/Enrichment/gene.txt")

# function
go_classify = function(gelst, type){
  goTable = data.frame(matrix(ncol = 9, nrow = 0)) %>%
    `colnames<-`(c("ID","Description","GeneRatio","BgRatio",
                   "pvalue","p.adjust","qvalue","geneID","Count"))
  
  pm1GO2Gene = filter(pm1GeneGo, Ontology == type) %>%
    select(GO, Gene)
  pm1GO2Term = filter(pm1GeneGo, Ontology == type) %>%
    select(GO, TERM)
  
  go_type = enricher(gelst,
                     pAdjustMethod = 'fdr',
                     minGSSize = 5,
                     TERM2GENE = pm1GO2Gene,
                     TERM2NAME = pm1GO2Term,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2)
  
  if (is.null(go_type)){
    go_type = goTable
  } else {
    go_type = data.frame(go_type@result,
                         stringsAsFactors = F)
  }
  return(go_type)
}

go_result = function(genelist){
  if (!exists("pm1GeneGo")) {
    message("Reading annotation...")
    pm1GeneGo = read_tsv("pm1.go")
  }
  geneKeep = genelist[genelist %in% pm1GeneGo$Gene]
  go_res = c("BP", "CC", "MF") %>% 
    map_df(~ go_classify(gelst = geneKeep, .x)) %>%
    as_tibble()
  return(go_res)
}
