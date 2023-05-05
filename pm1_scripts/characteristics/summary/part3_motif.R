#!/usr/bin/Rscript

# signature
## Author: yuki
## Date: 2023.04.08

## environment
rm(list = ls())
setwd("D://_CircRNA/results/10_revise/20230201/")
rltdir="D://_CircRNA/results/"

## packages
source("D://_Scripts/R/_LIBRARY.R")
library(paletteer)
library(rlang)
library(ggpubr)

# colors and themes
source("D://_Scripts/R/_COLOR/ggsci_color_theme.R")

# data
circ_info = read_tsv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circRNA_info_all.tsv"))
circ_flk_up = read_tsv(paste0(rltdir, "8_circ_WT_merge/6_flanking/3_motif/1_streme/flk.up/sequences.tsv")) %>%
    filter(`motif_P-value` < 0.05, grepl("pm1_", seq_ID)) %>%
    `colnames<-`(c("motif_struc", "motif_id", "p", "circ_id", "seq_score", "seq_class", "holdout")) %>%
    mutate(motif_type = "up")
circ_flk_down = read_tsv(paste0(rltdir, "8_circ_WT_merge/6_flanking/3_motif/1_streme/flk.down/sequences.tsv")) %>%
    filter(`motif_P-value` < 0.05, grepl("pm1_", seq_ID)) %>%
    `colnames<-`(c("motif_struc", "motif_id", "p", "circ_id", "seq_score", "seq_class", "holdout")) %>%
    mutate(motif_type = "down")

motif_bind = read_tsv(paste0(rltdir, "8_circ_WT_merge/6_flanking/3_motif/2_tomtom/motif_enrich.tsv"))

# body
circ_flk = rbind(circ_flk_up, circ_flk_down)
circ_flk %>% filter(motif_id == "STREME-1") %>% select(circ_id) %>% unique() %>% dim()

length(unique(circ_flk$circ_id))
2389/3891*100

table(circ_flk_up$circ_id %in% circ_flk_down$circ_id)
597/3891*100

# motif table



# hn and sf
motif_anno = data.frame(
    motif_struc = unique(c(hnRNP_motif, SF_motif)),
    motif_anno = c("hnRNP,SF","hnRNP","hnRNP,SF","hnRNP","hnRNP", rep("SF", 3))
) %>% left_join(motif_bind[grepl("HNRNP", motif_bind$gene_id),], by = c("motif_struc" = "motif")) %>%
    left_join(motif_bind[grepl("SF", motif_bind$gene_id),], by = c("motif_struc" = "motif"),
              suffix = c(".hnRNP", ".SF"))

hnRNP_motif = unique(unlist(motif_bind[grepl("HNRNP", motif_bind$gene_id),"motif"]))
SF_motif = unique(unlist(motif_bind[grepl("SF", motif_bind$gene_id),"motif"]))

circ_flk %>% select(circ_id) %>% unique() # 2389
hncirc = circ_flk %>% filter(motif_struc %in% hnRNP_motif) %>% select(circ_id) %>% unique() #852
sfcirc = circ_flk %>% filter(motif_struc %in% SF_motif) %>% select(circ_id) %>% unique() #1472
c(hncirc$circ_id, sfcirc$circ_id) %>% length()
intersect(hncirc$circ_id, sfcirc$circ_id) %>% length()

852/3891
835/3891
