#!/usr/bin/Rscript
setwd("D://_CircRNA/results/5_circRNA_merging/7_expression/")

library(tidyverse)
library(DESeq2)

# function
deFilter = function(res, id){
    type = c("res_na", "res_up", "res_down", "res_sig")
    idx_vars = paste(type, "idx", sep = "_")
    res_vars = paste(id, type, sep = "_")

    assign(idx_vars[1], is.na(res$pvalue) | is.na(res$padj) | is.na(res$log2FoldChange))
    assign(idx_vars[2], res$padj < 0.05 & res$log2FoldChange > 1)
    assign(idx_vars[3], res$padj < 0.05 & res$log2FoldChange < (-1))
    assign(idx_vars[4], res$padj < 0.05 & (abs(res$log2FoldChange) > 1))
    
    for (i in seq_along(type)){
        message(paste("Generating", res_vars[i], "by", idx_vars[i]), sep = " ")
        assign(res_vars[i], filter(res, get(idx_vars[i])), envir = .GlobalEnv)    
    }
}

deWrite = function(id){
    oldpath = getwd()
    dir.create(id)
    setwd(id)
    type = c("res_na", "res_up", "res_down", "res_sig")
    out_vars = paste(id, type, sep = "_")
    for (i in seq_along(type)){
        message(paste("Writing", out_vars[i], "to",  paste0(out_vars[i], ".tsv")))
        print(head(get(out_vars[i])))
        write_tsv(x = get(out_vars[i]), file = paste0(out_vars[i], ".tsv"))
    }
    setwd(oldpath)
}

############################### main ###############################
# read expression matrix

linear_ciriquant_trans = "_linear_by_CIRIquant/IN_merge_transcript_count_matrix.csv"
linear_ciriquant_gene = "_linear_by_CIRIquant/IN_merge_gene_count_matrix.csv"
linear_fcts_both_trans = "_linear_by_featureCounts/both/transcript_count_bothend.clean.tsv"
linear_fcts_both_gene = "_linear_by_featureCounts/both/gene_count_bothend.clean.tsv"

sfw = "fcts"

# read in
if (sfw == "ciriquant"){
  trans_count = read_csv(linear_ciriquant_trans) %>%
    as.data.frame() %>%
    column_to_rownames("transcript_id")
  
  gene_count = read_csv(linear_ciriquant_gene) %>%
    as.data.frame() %>%
    column_to_rownames("gene_id")
  
} else if (sfw == "fcts"){
  trans_count = read_tsv(linear_fcts_both_trans) %>%
    as.data.frame() %>%
    column_to_rownames("Geneid")
  
  gene_count = read_tsv(linear_fcts_both_gene) %>%
    as.data.frame() %>%
    column_to_rownames("Geneid")
} else {
  stop("Are you OK?")
}

message("*** Transcript Matrix ***")
head(trans_count)

message("*** Gene Matrix ***")
head(gene_count)

colData = data.frame("Sample" = colnames(gene_count), "Condition" = c(rep("Mycelium", 4), rep("Yeast", 4)))

message("*** ColData ***")
colData

# construct DESeq2 object
trans_dds = DESeqDataSetFromMatrix(countData = trans_count, colData = colData, design = ~ Condition)
trans_dds = DESeq(trans_dds)

gene_dds = DESeqDataSetFromMatrix(countData = gene_count, colData = colData, design = ~ Condition)
gene_dds = DESeq(gene_dds)

# differential expression results
de_trans_res = as.data.frame(results(trans_dds)) %>% 
    rownames_to_column("transcript_id")
de_gene_res = as.data.frame(results(gene_dds)) %>%
    rownames_to_column("gene_id")

deFilter(de_trans_res, "trans")
deFilter(de_gene_res, "gene")

deWrite("trans")
deWrite("gene")
