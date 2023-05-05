#!/usr/bin/Rscript

# signature
## Author: yuki
## Date: 2023.04.18

## environment
rm(list = ls())
dedir="D:/_CircRNA/results/9_WGCNA/1_RRHO/1_DEseq2/"
workdir="D://_CircRNA/results/10_revise/20230201/"
rltdir="D://_CircRNA/results/"
setwd(workdir)

## packages
source("D://_Scripts/R/_LIBRARY.R")
library(paletteer)
library(rlang)
library(ggpubr)

# colors and themes
source("D://_Scripts/R/_COLOR/ggsci_color_theme.R")

# function 

# data
circ_de = read_csv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circRNA_de.tsv")) %>%
    rename(circ_id_vertical = ...1)
circ_ratio_de = read_tsv(paste0(rltdir, "8_circ_WT_merge/7_expression/3_ratio/ratio/ratio_res.tsv")) %>%
    left_join(circ_info[,c(1,14)])
gene_de = read_tsv(paste0(rltdir, "8_circ_WT_merge/7_expression/1_linear/gene/gene_res.tsv"))

gene_exp = read_tsv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_gene_circ_CPG.csv"))
circ_bsj = read_csv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circRNA_bsj.csv"))

# body
# circ heatmap
de_circ_id = circ_de %>% filter(DE != 0)

ch_anno_col = data.frame(
    sample = factor(rep(c("mycelium", "yeast"), each = 3))
) %>% `rownames<-`(c(paste0("WTM",1:3), paste0("WTY",1:3)))

ch_ann_colors = list(
    sample = c(mycelium = col_nejm[2], yeast = col_nejm[1])
)

pdf(paste0(figdir, "Figure4/circ_heatmap_blue.pdf"), width = 3, height = 2.3)
gene_exp %>% filter(gene_id %in% de_circ_id$circ_id_vertical) %>%
    column_to_rownames("gene_id") %>%
    pheatmap::pheatmap(show_rownames = F,scale = "row",
                       #legend = FALSE, 
                       #annotation_legend = FALSE,
                       border=FALSE,
                       show_colnames = F,
                       annotation_col = ch_anno_col,
                       annotation_colors = ch_ann_colors,
                       colorRampPalette(c(col_nejm[2], "white",  col_nejm[1]))(50), 
                       treeheight_row = 15, treeheight_col = 15)
dev.off()

# go term
de_circ_go = read_tsv(paste0(rltdir, "8_circ_WT_merge/7_expression/4_enrichment/circ_de_go.tsv")) %>%
    arrange(pvalue) %>%
    mutate(list_hit = as.numeric(str_split(GeneRatio, pattern = "/", simplify = T)[,1]),
           list_bk = as.numeric(str_split(GeneRatio, pattern = "/", simplify = T)[,2]),
           FE = round(list_hit/list_bk, 2))
de_gene_go = read_tsv(paste0(rltdir, "8_circ_WT_merge/7_expression/4_enrichment/gene_de_go.tsv")) %>%
    arrange(pvalue) %>%
    mutate(list_hit = as.numeric(str_split(GeneRatio, pattern = "/", simplify = T)[,1]),
           list_bk = as.numeric(str_split(GeneRatio, pattern = "/", simplify = T)[,2]),
           FE = round(list_hit/list_bk, 4))


ggplot(de_circ_go, aes(x = FE,y = Description, color = p.adjust)) +
    geom_point() +
    scale_x_continuous(limits = c(0, 0.15), name = "fold enrichment") +
    scale_y_discrete(labels = \(x) {str_wrap(x, width = 25)}, name = "GO terms") +
    scale_color_continuous() +
    theme_bw() +
    mytheme +
    theme(
        axis.text.y = element_text(size = 7),
    )
ggsave(paste0(figdir, "Figure4/go_circ.pdf"), width = 80, height = 60, units = "mm")

de_gene_go %>% slice(1:6) %>%
ggplot(., aes(x = FE,y = Description, color = p.adjust)) +
    geom_point() +
    scale_x_continuous(limits = c(0, 0.12), name = "fold enrichment") +
    scale_y_discrete(labels = \(x) {str_wrap(x, width = 25)}, name = "GO terms") +
    scale_color_continuous() +
    theme_bw() +
    mytheme +
    theme(
        axis.text.y = element_text(size = 7)
    )
ggsave(paste0(figdir, "Figure4/go_gene.pdf"), width = 90, height = 60, units = "mm")

