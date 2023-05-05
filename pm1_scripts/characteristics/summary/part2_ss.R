#!/usr/bin/Rscript

# signature
## Author: yuki
## Date: 2023.04.07

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
circ_splice = read_tsv(paste0(rltdir, "8_circ_WT_merge/4_spliceSite/all.motifRev.txt"))
circ_marker = read_tsv(paste0(rltdir, "8_circ_WT_merge/5_structure/circ_marker_table_type.tsv"))

circ_ssite = read_tsv(paste0(rltdir, "8_circ_WT_merge/4_spliceSite/all.motifRev.table"), 
                      col_names = c("motif", "count", "freq")) %>%
    mutate(type = "circRNA")
gene_ssite = read_delim(paste0(rltdir, "8_circ_WT_merge/4_spliceSite/gene_intron.table"),
                      col_names = c("count", "freq", "motif"), delim = " ") %>%
    mutate(type = "gene",
           freq = -freq,
           count = -count) 

############################################################################################
# body
# length v1
circ_length_compare = circ_info %>% 
    select(circ_id_vertical, pred_len, len) %>%
    rename(D = len,
           P = pred_len) %>%
    mutate(length_type = ifelse(is.na(D), "IO", "NS"),
           len = ifelse(is.na(D), 0, D)) %>%
    pivot_longer(cols = P:D,
                 names_to = "len_type",
                 values_to = "length") %>%
    mutate(type = paste(length_type, len_type, sep = ""),
           type = factor(type, levels = c("IOP", "NSP", "NSD"))) %>%
    filter(length != 0) %>%
    left_join(circ_type, by = c("circ_id_vertical" = "circ_id"))

ggplot(circ_length_compare, aes(x = circ_type, y = log10(length+1), color = type)) +
    geom_boxplot(width  = 0.5, outlier.alpha = 0.5, outlier.size = 0.5, size = 0.5) +
    scale_color_manual(values = rev(col_nejm)[c(4,3,6)]) + 
    labs(x = "circ type") +
    theme_classic() +
    mytheme +
    theme(
        legend.position = "right",
        legend.background = element_rect(fill = NA),
        legend.title = element_blank()
    )
ggsave(paste0(figdir, "Figure3/lengthCompare.pdf"), width = 70, height = 48, units = "mm")


############################################################################################
# length v2
circ_info %>%
    mutate(type = ifelse(is.na(len), "Illumina", "Nanopore"),
           length =  ifelse(is.na(len), pred_len, len)) %>%
    ggplot(., aes(x = type, y = length, color = type, fill = type)) +
        geom_hline(yintercept = c(100, 1000, 10000),
                   linetype = "dashed", size = .15, color = "grey80") +
        geom_violin(size = 0.5, width = 0.5) +
        geom_boxplot(width = 0.2, fill = "white", size = 0.2,
                     outlier.color = NA) +
        scale_y_log10(name = "length (nt)") +
        scale_x_discrete(name = "") +
        scale_color_manual(values = col_nejm[5:6]) +
        scale_fill_manual(values = alpha(col_nejm[5:6]), .25) +
        theme_classic() +
        mytheme +
        theme(
            legend.position = "none"
        )

ggsave(paste0(figdir, "Figure3/lengthCompare_v2.pdf"), width = 55, height = 50, units = "mm")

circ_length_compare %>%
    filter(nano_supported == "yes") %>%
    select(circ_id_vertical, length, type) %>%
    pivot_wider(names_from = type, values_from = length) %>%
    ggplot(aes(x = NSD, y = NSP)) +
        geom_hline(yintercept = seq(0, 2000, 500),
               linetype = "dashed", size = .15, color = "grey80") +
        geom_point(color = col_nejm[6], size = 0.25) +
        geom_abline(slope = 1, size = 0.4, linetype = "dashed", color = "grey60") +
        scale_x_continuous(name = "length by Nanopore (nt)", expand = c(0, 0)) +
        scale_y_continuous(name = "length by Illumina (nt)", expand = c(0, 0)) +
        theme_classic() +
        mytheme
ggsave(paste0(figdir, "Figure3/lengthCompare_v3.pdf"), width = 55, height = 50, units = "mm")


# 1217 single
# 10378 - 1217 = 9161
# circularization of genes
circ_marker %>% filter(grepl("multi", circ_marker)) %>% select(gene_id) %>% unique() %>% dim() #1225
circ_marker %>% filter(grepl("single_exon_", circ_marker)) %>% select(gene_id) %>% unique() %>% dim() #47


data.frame(
    gene_type = c("single_exon", "multi_exons", "bg"),
    gene_number = c(1217, 9161, 0),
    circ_number = c(47, 1225, 0)
) %>% mutate("noncirc" = gene_number - circ_number) %>%
    pivot_longer(cols = circ_number:noncirc,
                 names_to = "type",
                 values_to = "number") %>%
    mutate(prop = number/gene_number *100,
           ypos = 0.75*prop,
           color = paste(gene_type, type, sep = "-")) %>%
    ggplot(., aes(x = gene_type, y = prop, fill = color)) +
        geom_bar(stat="identity", width=0.75, size = 0.25) +
        coord_polar("y", start=0)+
        geom_text(aes(y = ypos, label = paste0(round(prop, 2),"%")), color = "black", size=3) +
        scale_color_manual(values = col_jama) +
        scale_fill_manual(values = alpha(col_jama,  0.5)) +
        theme_void()

ggsave(paste0(figdir, "Figure3/SMgeneNumber.pdf"), width = 100, height = 60, units = "mm")

single_chi = matrix(c(1225, 47, 7936, 1170), byrow = T, nrow = 2) %>% chisq.test()
single_chi$p.value

# circRNA types
table(circ_marker$circ_marker)
# length diff in exons (intron retension)
circ_marker %>% filter(len_same == T, circ_marker == "multi_exon_exon") %>% 
    select(len, pred_len, compose, circ_marker)

data.frame(
    circ_type = c("seg", "se-meg", "si-meg", "me-meg", "ei-meg"),
    circ_number = c(52, 949, 34, 584, 191)
) %>% mutate(prop = circ_number/1810*1000,
             label = paste0(round(circ_number/1810*100, 2), "%")) %>%
    ggplot() +
        geom_hline(yintercept = seq(0, 1000, 250), size = 0.15, color = "grey80", linetype = "dashed") +
        geom_col(aes(x = circ_type, y = circ_number, fill = circ_type),
                 width = 0.25) +
        geom_text(aes(x = circ_type, y = circ_number, label = label), nudge_y = 75, size = 2) +
        scale_x_discrete(name = "") +
        scale_y_continuous(name = "circRNA number", limits = c(0, 1100), breaks = seq(0, 1100, 250)) +
        scale_fill_nejm() +
        scale_color_nejm() +
        theme_classic() +
        mytheme + 
        theme(
            axis.text.x = element_text(angle = 45),
            legend.position = "none"
        )
        
ggsave(paste0(figdir, "Figure3/typeNumber.pdf"), width = 75, height = 45, units = "mm")


# intron length
intron_number_length = read_tsv(paste0(rltdir, "8_circ_WT_merge/5_structure/intron_number_t01.bed"), 
                                col_names = c("chr", "ss", "es", "intron_id", "intron_length", "strand"))
circ_intron_feat = read_tsv(paste0(rltdir, "8_circ_WT_merge/5_structure/circ_feat_intron.tsv"), col_names = F)[,c(3,8)] %>%
    `colnames<-`(c("circ_id", "intron_id"))

gene_circ_intron = full_join(circ_intron_feat, intron_number_length) %>%
    mutate(feat_type = ifelse(is.na(circ_id), "other intron", "intron of circRNA"))

ggplot(gene_circ_intron, aes(x = log10(intron_length), color = feat_type)) +
    geom_density(size = 0.25) +
    scale_color_manual(values = col_nejm) +
    #scale_y_continuous(name = "") +
    scale_x_continuous(name = "log10(length of intron)") +
    theme_classic() +
    mytheme +
    theme(
        legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA,),
        panel.grid = element_blank()
    )
ggsave(paste0(figdir, "Figure3/intron_length1_clean.pdf"), width = 50, height = 45, units = "mm")

ggplot(gene_circ_intron, aes(x = feat_type, y = log10(intron_length), color = feat_type)) +
    geom_boxplot(width = 0.35, outlier.size = 0.15, outlier.alpha = 0.25, size = 0.25) +
    scale_color_manual(values = col_nejm) +
    scale_x_discrete(name = "") +
    scale_y_log10(name = "length (nt)") +
    theme_classic() +
    mytheme +
    theme(
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA)
    )

ggsave(paste0(figdir, "Figure3/intron_length2_clean.pdf"), width = 55, height = 40, units = "mm")

wilcox.test(filter(gene_circ_intron, feat_type == "circRNA")$intron_length,
            filter(gene_circ_intron, feat_type == "gene")$intron_length)

# splicing sites
circ_ssite = read_tsv(paste0(rltdir, "8_circ_WT_merge/4_spliceSite/all.motifRev.table"), 
                      col_names = c("motif", "count", "freq")) %>%
    mutate(label = ifelse(count > 4, motif, "others")) %>%
    group_by(label) %>%
    summarise(cnt = sum(count)) %>%
    mutate(type = "circRNA",
           freq = cnt/3891)
    
gene_ssite = read_delim(paste0(rltdir, "8_circ_WT_merge/4_spliceSite/gene_intron.table"),
                      col_names = c("cnt", "freq", "label"), delim = " ") %>%
    mutate(type = "gene",
           freq = -freq,
           cnt = -cnt,
           label = paste0(str_sub(label, 3, 4), str_sub(label, 1, 2))) %>%
    select(label, cnt, type, freq)

rbind(circ_ssite, gene_ssite) %>%
    mutate(motif = paste(str_sub(label, 3, 4), str_sub(label, 1, 2), sep = "-"),
           motif = ifelse(label == "others", "others", motif)) %>%
    filter(label != "AGGT") %>%
    ggplot(aes(x = motif, y = freq*100, fill = type)) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.1, color = "grey") +
    geom_hline(yintercept = seq(-3, 3, 1), linetype = "dashed", size = 0.1, color = "grey80") +
    geom_col(width = 0.3) +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "", limits = c(-3.5, 3.5),
                       breaks = seq(-3.5, 3.5, 1), 
                       labels = paste0(breaks = c(3.5, 2.5, 1.5, 0.5 ,seq(0.5, 3.5, 1)), "%")) +
    scale_fill_manual(values = col_jama[c(4,7)]) +
    theme_bw() +
    mytheme +
    theme(
        axis.text.x = element_text(angle = 45, size = 6),
        legend.position = c(0.8,0.25),
        legend.background = element_blank(),
        panel.grid = element_blank()
    )
ggsave(paste0(figdir, "Figure3/motif_non_others.pdf"), width = 55, height = 45, units = "mm")

rbind(circ_ssite, gene_ssite) %>%
    mutate(motif = paste(str_sub(label, 3, 4), str_sub(label, 1, 2), sep = "-"),
           motif = ifelse(label == "others", "others", motif)) %>%
    filter(label == "AGGT") %>%
    ggplot(aes(x = motif, y = freq*100, fill = type)) +
    geom_col(width = 0.5) +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "", breaks = seq(-100,100,40), 
                       labels = paste0(breaks = c(100,60,20,20,60,100), "%")) +
    scale_fill_manual(values = col_jama[c(4,7)]) +
    theme_bw() +
    mytheme +
    theme(
        axis.text.x = element_text(angle = 45, size = 6),
        legend.position = "none",
        legend.background = element_rect(fill = NA),
        panel.grid = element_blank()
    )
ggsave(paste0(figdir, "Figure3/motif_can_others.pdf"), width = 22, height = 45, units = "mm")

 
# match and dismatch
circ_marker %>% filter(grepl("multi_intron", circ_marker)) %>% select(circ_type, compose) %>% View()
multi_exact_table = circ_marker %>% 
    filter((grepl("multi_exon_exon", circ_marker))|(grepl("multi_intron", circ_marker))) %>% 
    filter(grepl("^m_*", compose)) %>% 
    filter(grepl("*_m$", compose))

circ_marker_match = circ_marker %>%
    mutate(match_type = ifelse(grepl("match", circ_marker), "match", "dismatch"),
           match_type = ifelse(circ_id_vertical %in% multi_exact_table$circ_id_vertical, 
                          "match", match_type)
           ) %>%
    left_join(circ_splice[,c(5,10)], by = c("circ_id_vertical" = "circID")) %>%
    mutate(
        motif_type = ifelse(motifRev == "AGGT", "canonical", "non-canonical")
    )

# canonical and non-canonical vs match
table(circ_marker_match$match_type, circ_marker_match$motif_type) %>% chisq.test()

# non-canonical site plot
names(table(circ_marker_match$motifRev)[table(circ_marker_match$motifRev) > 5])
circ_marker_match_non = circ_marker_match %>%
    filter(motifRev %in% names(table(circ_marker_match$motifRev)[table(circ_marker_match$motifRev) > 5])) %>%
    filter(motif_type == "non-canonical") %>%
    select(gene_id, compose, circ_marker, match_type, motifRev)

circ_marker_match_non %>% filter(motifRev == "AGAT", match_type == "match")

circ_marker_match_non %>% 
    select(match_type, motifRev) %>%
    table() %>%
    as.data.frame() %>%
    mutate(label = paste(str_sub(motifRev, 3, 4), str_sub(motifRev, 1, 2), sep = "-")) %>%
    ggplot() +
        geom_hline(yintercept = seq(10, 30, 10), linetype = "dashed", size = 0.1, color = "grey80") +
        geom_col(aes(x = label, y = Freq, fill = match_type),
                 width = 0.35) +
        scale_x_discrete(name = "non canonical motif") +
        scale_y_continuous(name = "number", expand = c(0, 0)) +
        scale_fill_manual(values = col_jama[5:6]) +
        theme_classic() +
        mytheme +
        coord_flip()

ggsave(paste0(figdir, "Figure3/motif_match.pdf"), width = 75, height = 45, units = "mm")

# exon and intron circRNA in dismatch
match_chi = circ_marker_match %>%
    filter(motif_type == "non-canonical") %>%
    mutate(intron_type = ifelse(grepl("intron", circ_marker), "intron", "exon")) %>%
    select(match_type, intron_type) %>%
    table() %>%
    chisq.test()
match_chi$p.value

