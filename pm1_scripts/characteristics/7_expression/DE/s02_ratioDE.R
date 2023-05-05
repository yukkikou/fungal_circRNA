#!/usr/bin/Rscript

setwd("D://_CircRNA/results/5_circRNA_merging/7_expression/ratio/")

library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(apeglm)
library(fdrtool)

#rm(list = ls())

# function
matrixPrepare = function(circ_bsj, circ_info, gene_cnt){
  # only focus on genes product single circRNA
  multiIdx = base::grep(circ_info$gene_id, pattern = ",")
  circSingle = circ_info$gene_id[-multiIdx]
  circSingle = circSingle[!is.na(circSingle)]
  circGeneFre = table(circSingle)
  
  circ_info_single = circ_info %>% filter(gene_id %in% circSingle) %>%
    select(circ_id, gene_id)
  circ_bsj_single = circ_bsj %>% right_join(circ_info_single) %>%
    as.data.frame() %>%
    arrange(gene_id) %>%
    mutate(gene_name = NA)
  for (i in 1:length(circGeneFre)){
    if (circGeneFre[i] == 1){
      gene_tmp = str_split(names(circGeneFre[i]), pattern = "\\.")[[1]][1]
      circ_bsj_single[circ_bsj_single$gene_id == names(circGeneFre[i]), "gene_name"] = paste0(gene_tmp, ".c1")
    } else {
      gene_tmp = str_split(names(circGeneFre[i]), pattern = "\\.")[[1]][1]
      gene_num = circGeneFre[i]
      for (j in 1:gene_num){
        circ_bsj_single[circ_bsj_single$gene_id == names(circGeneFre[i]), "gene_name"][j] = paste0(gene_tmp, ".c", j)
      }
    }
  }
  gene_cnt_single = gene_cnt %>% right_join(select(circ_bsj_single, gene_id, gene_name),
                                            by = c("Geneid" = "gene_id"))
  message("Generate circ_bsj_matrix...")
  circ_bsj_matrix <<- circ_bsj_single %>% 
    select(IN_M1:IN_Y4,gene_name) %>%
    column_to_rownames("gene_name")
  print(head(circ_bsj_matrix))
  
  message("Generate gene_cnt_matrix...")
  gene_cnt_matrix <<- gene_cnt_single %>%
    select(C1:gene_name) %>%
    column_to_rownames("gene_name") %>%
    `colnames<-`(c(paste0("IN_M",1:4), paste0("IN_Y",c(3:4,1:2))))
  print(head(gene_cnt_matrix))
  
  message("Generate circ_to_gene...")
  circ_to_gene <<- circ_bsj_single %>% select(circ_id, gene_id, gene_name)
  print(head(circ_to_gene))
}

DESeq2Rex <- function(rnaCntTable, riboCntTable, rnaCond, riboCond,
                      contrast=NULL, minMeanCount=0.1) {
  
  ### append sample name with .rna or .ribo
  ### in case user uses sample sample name
  colnames(rnaCntTable) <- paste0(colnames(rnaCntTable), ".rna")
  colnames(riboCntTable) <- paste0(colnames(riboCntTable), ".ribo")
  
  ### input validation
  if (!identical(rownames(rnaCntTable), rownames(riboCntTable)))
    stop ("RNA- and Ribo-seq data must have the same set of genes")
  
  if (!is.data.frame(rnaCond)) rnaCond <- data.frame(cond = rnaCond)
  if (!is.data.frame(riboCond)) riboCond <- data.frame(cond = riboCond)
  
  if (!identical(colnames(rnaCond), colnames(riboCond)))
    stop("RNA- and Ribo-seq data must have the same set of conditions")
  
  if (ncol(rnaCntTable) != nrow(rnaCond))
    stop(paste("RNA-seq count table must have the",
               "same number of samples as in rnaCond"))
  
  if (ncol(riboCntTable) != nrow(riboCond))
    stop(paste("Ribo-seq count table must have the",
               "same number of samples as in riboCond"))
  
  # if (minMeanCount < 1)
  #   stop("minMeanCount must at least be 1")
  
  ### filter out low read count
  keep.rna <- rownames(rnaCntTable)[rowMeans(rnaCntTable) >= minMeanCount]
  keep.ribo <- rownames(riboCntTable)[rowMeans(riboCntTable) >= minMeanCount]
  keep <- intersect(keep.rna, keep.ribo)
  rnaCntTable <- rnaCntTable[keep,]
  riboCntTable <- riboCntTable[keep,]
 
  message(paste("Remaining circRNAs are", dim(rnaCntTable)[1]))
   
  numCond <- ncol(rnaCond)
  numRNASmps <- nrow(rnaCond)
  numRiboSmps <- nrow(riboCond)
  
  ### combine counts
  combCntTbl <- cbind(rnaCntTable, riboCntTable)
  
  message("combining design matrix")
  
  combinedCond <- rbind(rnaCond, riboCond)
  combinedCond <- combinedCond[,rep(1:ncol(combinedCond),2)]
  INTERCEPT <- c(rep("CONTROL", numRNASmps), rep("TREATED", numRiboSmps))
  combinedCond <- cbind(combinedCond[1:numCond], INTERCEPT,
                        combinedCond[(numCond+1):ncol(combinedCond)])
  for( i in (numCond+2) : ncol(combinedCond)) {
    combinedCond[1:numRNASmps,i] <- combinedCond[1,i]
  }
  colnames(combinedCond)[(numCond+2):ncol(combinedCond)] <- paste0("EXTRA",
                                                                   seq(numCond))
  extendedConds <- colnames(combinedCond)
  fmla <- as.formula(paste("~", paste(extendedConds, collapse= "+")))
  dds <- DESeqDataSetFromMatrix(countData = combCntTbl,
                                colData = combinedCond,
                                design = fmla)
  
  message("applying DESeq2 to modified design matrix")
  
  ## apply new design matrix with combined count table to DESeq2
  dds <- DESeq(dds)
  if(is.null(contrast)) {
    res <- results(dds)
  } else {
    contrast[1] <- paste0("EXTRA", which(colnames(combinedCond)==contrast[1]))
    res <- results(dds, contrast=contrast)
  }
  
  ## order results by gene names
  res <- res[order(res$padj),]
  print(head(res))
  res
}

## main
# read in

circ_bsj_file = "IN_merge_circRNA_bsj.csv"
circ_de_file = "IN_merge_circRNA_de.tsv"
circ_info_file = "IN_merge_circRNA_info.csv"
gene_cnt_file = "_linear_by_featureCounts/both/gene_count_bothend.clean.tsv"

circ_bsj = read_csv(circ_bsj_file)
circ_info = read_csv(circ_info_file)
gene_cnt = read_tsv(gene_cnt_file)

de_circ_res = read_csv(circ_de_file)
colnames(de_circ_res)[1] = "circ_id"

# generate dds
matrixPrepare(circ_bsj, circ_info, gene_cnt)
de_ratio_res = DESeq2Rex(rnaCntTable = gene_cnt_matrix,
          riboCntTable = circ_bsj_matrix,
          rnaCond = c(rep("Mycelium",4),rep("Yeast",4)),
          riboCond = c(rep("Mycelium",4),rep("Yeast",4))) %>% 
  as.data.frame() %>% rownames_to_column("circtrans_id")

# output
deFilter(de_ratio_res, "ratio")
deWrite("ratio")
write_tsv(circ_to_gene, "circ_to_gene.tsv")


# special check for circ-DS1
de_ratio_res[grep(de_ratio_res$circtrans_id, pattern = "TM020485"),]
circ_to_gene %>% filter(gene_id == "TM020485.00")

# visualization
de_ratio_res %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  geom_point(aes(size=-log10(padj), color= -log10(padj)))+
  labs(x = "log2(FoldChange)", y = "-log10(FDR)") +
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  scale_size_continuous(range = c(1,3))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=10),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        legend.text = element_text(size=6,color = "black"),
        legend.title = element_text(size=8),
        legend.key.size = unit(11, "pt")) +
  theme(panel.grid = element_blank())
