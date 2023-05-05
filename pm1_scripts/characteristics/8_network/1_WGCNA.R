#!/usr/bin/Rscript

# Signature
## Author: yuki
## Date: 2023.02.09

# packages
library(WGCNA)
library(reshape2)
library(tidyverse)
library(data.table)
library(ggsci)

# environment
rm(list = ls())
rltdir = "D://_CircRNA/results/7_mutants/"
workdir = "D://_CircRNA/results/9_WGCNA/2_WGCNA/"
setwd(workdir)

enableWGCNAThreads()
ALLOW_WGCNA_THREADS=4
#disableWGCNAThreads()
#memory.limit(size = 20000)

# colors and theme
source("D://_CircRNA/scripts/color_theme.R")

# source
source("D://_CircRNA/scripts/Enrichment/s00_goFunction.R")

# function
source("D://_CircRNA/scripts/network/WGCNA_function.R")

wgcnaOnestep = function(con, gcCpmFV){
  datExpr0 = as.data.frame(t(gcCpmFV[,-1])) %>%
    `colnames<-`(gcCpmFV$feat_id)
  
  datExpr = datExprClean(expr = datExpr0)
  
  # samples cluster
  col_index = get(paste("col",  strsplit(con, "")[[1]][1], "index", sep = "_"))
  
  sampleTree(expr = datExpr,
             traitColors = rep(col_nejm[1:3], each = 3), 
             anno = colData$phase[col_index])
  
  # select soft threshold power
  sft = powerSelect(expr = datExpr)
  
  # dynamic motif construction
  if (is.na(sft)){
    sft = 20  
  } else {
    message(paste0("Using softthred ", sft))
  }
  
  lst = dyMotif(expr = datExpr, softPower = sft)
  
  # connect to colData and keep significant modules
  mtCor = toColData(expr = datExpr, mcol = lst$moduleColors, ns = 9, anno = colDataMtx[col_index,])
  write_tsv(mtCor, paste0(con, "_module_train_cor_pvalue.tsv"))
  
  mtCor_sig = mtCor %>%
    filter(abs(phase.cor) > 0.5, abs(phase.pvalue) < 0.2) %>%
    as.data.frame() %>%
    .[,"module"] %>%
    str_split(., pattern = "ME", simplify = T) %>%
    .[,2]
  
  print(mtCor_sig)
  
  # gene and module
  geneMod = data.frame(
    "feat_id" = colnames(datExpr),
    "module" = paste(con, moduleColors, sep = "_")
  )
  
  write_tsv(geneMod, paste(con, "geneMod.tsv", sep = "_"))
  
  # gene in interested module
  for (i in seq_along(mtCor_sig)){
    message(paste0("Now is module ", mtCor_sig[i]))
    #calGeneModuleColdata(anno = colDataMtx$phase[col_m_index], prefix = "phase", moduleCol = mtCor_sig[i])
    
    exGeneFromModule(expr = datExpr, mod_cols = lst$moduleColors, moduleCol = mtCor_sig[i]) %>%
      write_lines(paste(con, mtCor_sig[i], "modules_gene.tsv", sep = "_"))
    
    # output module
    outModuleToCyto(expr = datExpr, tom = lst$TOM, mod_cols = lst$moduleColors,
                    moduleCol = mtCor_sig[i], filter = F)
    
    # go enrichment
    # yeast
    tmp_go = exGeneFromModule(expr = datExpr, mod_cols = lst$moduleColors, 
                              moduleCol = mtCor_sig[i]) %>%
      go_result() %>%
      filter(p.adjust < 0.05)
    
    write_tsv(tmp_go, paste(con, mtCor_sig[i], "modules_enrich.tsv", sep = "_"))
    assign(paste(con, mtCor_sig[i], "go_enrich", sep = "_"), tmp_go, envir = .GlobalEnv)
  }
  circDS_col = exGeneModuleCol(expr = datExpr, mod_cols = lst$moduleColors, 
                  feat = "pm1_02:1328636|1328808", search = F)
  message(paste("CircDS module is", circDS_col))
  outModuleToCyto(expr = datExpr, tom = lst$TOM, mod_cols = lst$moduleColors,
                  moduleCol = circDS_col, filter = F)
  
}

## read in
# sample information
colData = data.frame(
  mutant = c(rep("KD", 6), rep("OE", 6), rep("WT", 6)),
  condition = rep(rep(c("mycelium", "yeast"), each = 3), 3),
  sample = c(paste0("KD", 2:4, "M"),paste0("KD", 2:4, "Y"),
             paste0("OE", LETTERS[4:6], "M"), paste0("OE", LETTERS[4:6], "Y"),
             paste0("WTM", 1:3),paste0("WTY", 1:3)),
  phase = rep(1:6, each = 3)
) %>%
  mutate(label = paste(mutant, condition, sep = "-"))

colDataMtx = colData %>% 
  column_to_rownames("sample") %>%
  mutate(
    sample = colData$sample
  )

for (i in 1:dim(colDataMtx)[2]){
  colDataMtx[,i] = as.numeric(factor(colDataMtx[,i]))
}
colDataMtx$sample = 1:18

# expression counts matrix
geneMtx = read_tsv(paste0(rltdir, "2_linear/linear_gene_count.clean.tsv"))
col_geneMtx = str_split(colnames(geneMtx)[-1], pattern = "\\/", simplify = T)[,2]
colnames(geneMtx) = c("feat_id", col_geneMtx)
geneMtx = geneMtx[, c("feat_id", colData$sample)]

circMtx = read_csv(paste0(rltdir, "1_circ/merge_circRNA_bsj.csv")) %>%
  rename(feat_id = circ_id)
circMtx = circMtx[, c("feat_id", colData$sample)]
#save(list = ls(), file = "exp_colData_raw.RData")

## merge expression matrix (all)
gcMtx = rbind(geneMtx, circMtx)

# clean expression matrix
# mycelium and yeast FV expression matrix separately
exp_m_index = c(1:4,8:10,14:16)
exp_y_index = c(1,5:7,11:13,17:19)
col_m_index = c(1:3,7:9,13:15)
col_y_index = c(4:6,10:12,16:18)

for (i in c("gene", "circ")){
  for (j in c("m", "y")){
    mtx_tmp = get(paste0(i, "Mtx"))
    e_idx = get(paste("exp", j, "index", sep = "_"))
    c_idx = get(paste("col", j, "index", sep = "_"))
    tmp = top1000zs(mtx = mtx_tmp[,e_idx], thd = 1, 
              colData = colData[c_idx,], 
              prec = 1, number = 3000)
    assign(paste0(i, "CpmFV1000", toupper(j)), tmp)
  }
}

gcCpmFVM = rbind(circCpmFV1000M, geneCpmFV1000M) %>%
  as.data.frame() %>%
  rownames_to_column("feat_id") %>%
  filter(feat_id %in% c(rownames(circCpmFV1000M), rownames(geneCpmFV1000M)))

gcCpmFVY = rbind(circCpmFV1000Y, geneCpmFV1000Y) %>%
  as.data.frame() %>%
  rownames_to_column("feat_id") %>%
  filter(feat_id %in% c(rownames(circCpmFV1000Y), rownames(geneCpmFV1000Y)))

# transfor to matrix and clean
getwd()
dir.create("mycelium")
setwd("mycelium")
wgcnaOnestep(con = "mycelium", gcCpmFV = gcCpmFVM)
setwd("../")

dir.create("yeast")
setwd("yeast")
wgcnaOnestep(con = "yeast", gcCpmFV = gcCpmFVY)
setwd("../")

for (i in grep("go_enrich", ls(), value = T)){
  tmp = get(i)
  print("****************")
  print(i)
  print(tmp$Description)
}

# circDS-1 coexpression
annoCyto = function(file){
    circ_id = "pm1_02:1328636|1328808"
    circ_name = "circDS-1"
    read_tsv(file) %>%
        select(fromNode:direction) %>%
        mutate(from_keep = grepl(fromNode, pattern = circ_id),
               to_keep = grepl(toNode, pattern = circ_id),
               final_keep = from_keep | to_keep) %>%
        filter(final_keep) %>%
        arrange(-weight) %>%
        left_join(pm1ReportedGene, by = c("fromNode" = "gene_id")) %>%
        left_join(pm1ReportedGene, by = c("toNode" = "gene_id"),
                  suffix = c(".from", ".to")) %>%
        left_join(pm1Omics, by = c("fromNode" = "SequenceName")) %>%
        left_join(pm1Omics, by = c("toNode" = "SequenceName"),
                  suffix = c(".from", ".to")) %>%
        rename(fromAltName = gene_name.from,
               toAltName = gene_name.to) %>%
        mutate(fromAltName = ifelse(fromNode == circ_id, circ_name, fromAltName),
               toAltName = ifelse(toNode == circ_id, circ_name, toAltName))
}


mycelium_circDS_edge = annoCyto("mycelium/CytoscapeInput-edges-grey60_circDS.txt") %>% head(30) %>% 
    mutate(toAltName = SequenceDescription.to,
           toAltName = ifelse(is.na(toAltName), toNode, toAltName)) %>%
    select(fromNode:direction, fromAltName:toAltName)

write_tsv(mycelium_circDS_edge, "mycelium/Final_mycelium_circDS_edge_30.tsv")

yeast_circDS_edge = annoCyto("yeast/CytoscapeInput-edges-lightgreen_circDS.txt") %>% head(30) %>% 
    mutate(toAltName = SequenceDescription.to,
           toAltName = ifelse(is.na(toAltName), toNode, toAltName)) %>%
    select(fromNode:direction, fromAltName:toAltName)
write_tsv(yeast_circDS_edge, "yeast/Final_yeast_circDS_edge_30.tsv")


# output
write_csv(annoCyto("mycelium/CytoscapeInput-edges-grey60_circDS.txt"),
          "mycelium/mycelium_circDS_edge.tsv")
write_csv(annoCyto("yeast/CytoscapeInput-edges-lightgreen_circDS.txt"),
          "yeast/yeast_circDS_edge.tsv")
