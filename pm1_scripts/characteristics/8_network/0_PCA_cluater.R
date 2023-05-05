#!/usr/bin/Rscript

## environment
rm(list = ls())
setwd("D://_CircRNA/results/9_WGCNA/")
rltdir="D://_CircRNA/results/7_mutants/"

## packages
library(tidyverse)
library(ggsci)
library(patchwork)

# colors
show_col(pal_nejm("default")(8))
show_col(pal_jama("default")(7))
col_nejm = c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF",
             "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF")
col_jama = c("#374E55FF", "#DF8F44FF", "#00A1D5FF", "#B24745FF",
             "#79AF97FF", "#6A6599FF", "#80796BFF")

# themes
mytheme = theme(plot.title = element_text(hjust = 0.5),
                text=element_text(size=10),
                axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
                axis.text.y = element_text(size = 8,color = "black"),
                legend.text = element_text(size=6,color = "black"),
                legend.title = element_text(size=8),
                legend.key.size = unit(11, "pt"))

## function
callog2Cpm = function(mtx){
  mtx_tmp = mtx[, -1]
  libSize = colSums(mtx_tmp)
  mtx_out = data.frame(matrix(nrow = dim(mtx)[1], ncol = dim(mtx)[2]))
  for (i in 1:dim(mtx_tmp)[2]){
    mtx_out[, (i+1)] = log2(mtx_tmp[, i]/libSize[i]*1000000+1)
  }
  print(colSums(mtx_out))
  mtx_out[, 1] = mtx[, 1]
  colnames(mtx_out) = colnames(mtx)
  as_tibble(mtx_out)
}

colFpkm = function(){
  
}

mtxFilter = function(mtx, thd, colData, prec){
  mtx_tmp = mtx %>%
    pivot_longer(cols = colnames(mtx)[-1], 
                 names_to = "sample",
                 values_to = "count") %>%
    left_join(colData, by = c("sample" = "sample")) %>%
    mutate(keep = as.numeric(count > thd)) %>%
    group_by(feat_id, phase) %>%
    summarise(k_score = sum(keep)/3,
              k_keep = as.numeric(k_score >= prec)) %>%
    group_by(feat_id) %>%
    summarise(f_score = sum(k_keep),
              f_keep = (f_score) > 0) %>%
    filter(f_keep == T) %>%
    as.data.frame()
  
  mtx %>% filter(feat_id %in% mtx_tmp$feat_id)
}

top1000zs = function(mtx, thd, colData, prec, number){
  mtx_tmp = mtx
  thd_tmp = thd
  colData_tmp = colData
  prec_tmp = prec
  # expression normalization  
  cpm = callog2Cpm(mtx)
  
  # expression filtering and low variance
  cpmFvars = mtxFilter(mtx = mtx_tmp, thd = thd_tmp, colData = colData_tmp, prec = prec_tmp) %>%
    column_to_rownames("feat_id") %>%
    apply(., 1, var)
  
  keep_feat = names(sort(cpmFvars, decreasing = T)[1:number])
  
  # z-score
  keep_cpm = cpm %>% as.data.frame() %>%
    filter(feat_id %in% keep_feat) 
  
  keep_cpm %>%
    column_to_rownames("feat_id") %>%
    t() %>%
    apply(., 1, function(x) scale(x , center = T, scale = T)) %>%
    `rownames<-`(keep_cpm$feat_id)
}

pcaPlot = function(id_mtx, prefix, colData, leading_number = 100){
  PC = prcomp(t(id_mtx))
  message(paste0("Generate variable: ", prefix, "_pca"))
  assign(paste0(prefix, "_pca"), PC, envir = .GlobalEnv)
  
  PC_sd <- setNames(PC$sdev , paste0("PC",1:length(PC$sdev)))
  PC_var_expl <- round(( PC_sd^2 ) / sum(PC_sd^2) * 100, 2)
  
  xlab = paste0("PC1 ","(", paste0(as.character(PC_var_expl[1]), "%", ")", sep=""))
  ylab = paste0("PC2 ","(", paste0(as.character(PC_var_expl[2]), "%", ")", sep=""))
  if (grepl(prefix, pattern = "circ")){
    cols = col_nejm[1:3]
  } else {
    cols = col_nejm[4:6]
  }
  
  plot_pca = PC$x[,1:2] %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    left_join(colData) %>%
    ggplot(aes(x = PC1, y = PC2, color = label, label = sample)) +
    geom_point(size = 1.5) +
    geom_text(nudge_x = 1, nudge_y = 1, check_overlap = T, size = 2.5) +
    labs(x = xlab, y = ylab) +
    scale_color_manual(values = cols) +
    theme_bw() +
    mytheme +
    theme(
      legend.position = "top",
      legend.title = element_blank()
    )
  
  message(paste0("Generate variable: ", "plot_", prefix, "_pca"))
  assign(paste0("plot_", prefix, "_pca"), plot_pca, envir = .GlobalEnv)
  
  leading_genes = PC$rotation
  
  leading_PC1 = sort(leading_genes[,1],decreasing = T)[1:leading_number] %>%
    as.data.frame() %>%
    rownames_to_column("feat_id") %>%
    `colnames<-`(c("feat_id", "rotation")) %>%
    mutate(pc = "PC1")
  leading_PC2 = sort(leading_genes[,2],decreasing = T)[1:leading_number] %>%
    as.data.frame() %>%
    rownames_to_column("feat_id") %>%
    `colnames<-`(c("feat_id", "rotation")) %>%
    mutate(pc = "PC2")
  
  leading_PC = rbind(leading_PC1, leading_PC2) %>% as_tibble()
  message(paste0("Generate variable: ", prefix, "_pca_leadingGene"))
  assign(paste0(prefix, "_pca_leadingGene"), leading_PC, envir = .GlobalEnv)
}

# example
# geneCpmFV1000M = top1000zs(mtx = geneMtx[,exp_m_index], thd = 1, 
#                            colData = colData[col_m_index,], 
#                            prec = 1, number = 1000)
# pcaPlot(id_mtx = geneCpmFV1000M, prefix = "geneM", 
#         colData = colData[col_m_index,], leading_number = 20)

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

# expression counts matrix
geneMtx = read_tsv(paste0(rltdir, "2_linear/linear_gene_count.clean.tsv"))
col_geneMtx = str_split(colnames(geneMtx)[-1], pattern = "\\/", simplify = T)[,2]
colnames(geneMtx) = c("feat_id", col_geneMtx)
geneMtx = geneMtx[, c("feat_id", colData$sample)]

circMtx = read_csv(paste0(rltdir, "1_circ/merge_circRNA_bsj.csv")) %>%
  rename(feat_id = circ_id)
circMtx = circMtx[, c("feat_id", colData$sample)]

circRatio = read_csv(paste0(rltdir, "1_circ/merge_circRNA_ratio.csv")) %>%
  rename(feat_id = circ_id)
circRatio = circRatio[, c("feat_id", colData$sample)]

# expression merge
gcMtx = rbind(geneMtx, circMtx)

## body
# PCA
exp_m_index = c(1:4,8:10,14:16)
exp_y_index = c(1,5:7,11:13,17:19)
col_m_index = c(1:3,7:9,13:15)
col_y_index = c(4:6,10:12,16:18)

for (i in c("gene", "circ")){
  for (j in c("m", "y")){
    mtx_tmp = get(paste0(i, "Mtx"))
    e_idx = get(paste("exp", j, "index", sep = "_"))
    c_idx = get(paste("col", j, "index", sep = "_"))
    top1000zs(mtx = mtx_tmp[,e_idx], thd = 1, 
                               colData = colData[c_idx,], 
                               prec = 1, number = 1000) %>%
      pcaPlot(id_mtx = ., prefix = paste0(i,j), 
              colData = colData[c_idx,], leading_number = 20)
  }
}

(plot_circm_pca + plot_circy_pca) / (plot_genem_pca + plot_geney_pca)
ggsave("pca_circ_gene.pdf", height = 120, width = 120, units = "mm")


# heatmap
geneCpmFV1000 = top1000zs(mtx = geneMtx, thd = 1, colData = colData, prec = 1, number = 1000)
circCpmFV1000 = top1000zs(mtx = circMtx, thd = 1, colData = colData, prec = 1, number = 1000)
#circCpmFV1000[(rownames(circCpmFV1000) == "pm1_02:1328636|1328808"),]

anno_col = list(
  condition = c(yeast = col_nejm[1], mycelium = col_nejm[2]),
  mutant = c(KD = col_nejm[3], OE = col_nejm[4], WT = col_nejm[5])
)
pheatmap(geneCpmFV1000, scale = "none", 
         show_rownames = F, 
         clustering_method = "median", 
         treeheight_row = 10, 
         treeheight_col = 10, 
         annotation_col = column_to_rownames(colData[,1:3], "sample"),
         annotation_colors = anno_col)

pheatmap(circCpmFV1000, scale = "none", 
         show_rownames = F, 
         clustering_method = "median", 
         treeheight_row = 30, 
         treeheight_col = 10, 
         annotation_col = column_to_rownames(colData[,1:3], "sample"),
         annotation_colors = anno_col)

# cor and distance
# cor only
cor(geneCpmFV1000, method = "spearman") %>%
  pheatmap(., scale = "none")
cor(circCpmFV1000, method = "spearman") %>%
  pheatmap(., scale = "none")

# cor to distance
cor_distance = -(cor(geneCpmFV1000) -1)/2
cor_distance %>%  
  as.dist() %>%
  pheatmap(., scale = "none")

cor_distance = -(cor(circCpmFV1000) -1)/2
cor_distance %>%  
  as.dist() %>%
  pheatmap(., scale = "none")

# distance only
geneCpmFV1000 %>% t() %>%
  as.data.frame() %>%
  `rownames<-`(colnames(geneCpmFV1000)) %>%
  dist(., method = "manhattan") %>%
  pheatmap()

circCpmFV1000 %>% t() %>%
  as.data.frame() %>%
  `rownames<-`(colnames(geneCpmFV1000)) %>%
  dist(., method = "manhattan") %>%
  pheatmap()

## output
