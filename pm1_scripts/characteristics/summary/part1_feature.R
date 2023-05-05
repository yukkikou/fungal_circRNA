##!/usr/bin/Rscript

# signature
## Author: yuki
## Date: 2023.04.06

## environment
#rm(list = ls())
workdir="D://_CircRNA/results/10_revise/20230201/"
rltdir="D://_CircRNA/results/"
figdir = "D://_CircRNA//figure//_panel//revise//20230210/"
setwd(workdir)

## packages
source("D://_Scripts/R/_LIBRARY.R")
library(paletteer)
library(rlang)
library(ggpubr)
#install.packages("ggExtra")
#install.packages("aplot")
library(ggExtra)

# colors and themes
source("D://_Scripts/R/_COLOR/ggsci_color_theme.R")

# functions

# data
# for number and isoform compare
circ_info = read_tsv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circRNA_info_all.tsv"))
human_atlas_gene = read_delim("D://_CircRNA/database/circAtlas_human/human_circ.lst", 
                              col_names = c("bl","circ_number", "gene_id"), delim = " ")[,2:3]
circ_gene_gtf = read_delim(paste0(rltdir, "8_circ_WT_merge/merge/circ_gene_single.lst"), 
                         col_names = c("bl","circ_number", "gene_id"), delim = " ")[,2:3]
linear_gene_gtf = read_tsv("D://_IGV/Talaromyces_marneffei/Annotation/gene_linear.number",
                           col_names = c("linear_number","gene_id"))
#inner_join(circ_gene_gtf, linear_gene_gtf)

# for expression
circ_bsj = read_csv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circRNA_bsj.csv"))
circ_bsj_long = read_csv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circRNA_bsj_longer.csv"))
circ_ratio = read_csv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circRNA_ratio.csv"))
circ_ratio_long = read_csv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circRNA_ratio_longer.csv"))
circ_cpg = read_tsv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circ_CPG.tsv"))
circ_cpg_long = read_tsv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circ_CPG_longer.tsv"))

gene_exp = read_tsv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_gene_circ_CPG.csv"))

########################################################################################################
# body
# isoform number compared with human
quantile(circ_gene_gtf$circ_number, probs = seq(0, 1 ,0.05))
mean(circ_gene_gtf$circ_number)
quantile(human_atlas_gene$circ_number, probs = seq(0, 1 ,0.05))

# isoform plot
data.frame(
    "type" = c("genes deriving circRNA", "genes not deriving circRNA"),
    "number" = c(2025, 8353)
) %>%
    mutate(type = factor(type, levels = c("genes not deriving circRNA", "genes deriving circRNA"))) %>%
    mutate(prop = number / sum(number) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop ) %>%
    ggplot(., aes(x = "", y = prop, color = type, fill = type)) +
        geom_bar(stat="identity", width=1, size = 0.25) +
        coord_polar("y", start=0)+
        geom_text(aes(y = ypos, label = paste0(round(prop, 2),"%")), color = "black", size=3) +
        scale_color_manual(values = col_nejm[c(4,4)]) +
        scale_fill_manual(values = c(alpha(col_nejm[4], 0.2), alpha(col_nejm[4], 0.7))) +
        theme_void() +
        theme(
            legend.position = "bottom",
            legend.title = element_blank()
        )
ggsave(paste0(figdir, "Figure2/pie.pdf"), width = 80, height = 60, units = "mm")

table(circ_gene_gtf$circ_number) %>% as.data.frame() %>%
    `colnames<-`(c("isonumber", "freq")) %>%
    mutate(isonumber = as.numeric(isonumber),
           grp = ifelse(isonumber < 10, isonumber, ">9")) %>%
    group_by(grp) %>%
    summarise(number = sum(freq)) %>%
    mutate(grp = factor(grp, levels = c(1:9, ">9"))) %>%
    ggplot(., aes(x = grp, y = number)) +
    geom_col(width = 0.4, fill = col_nejm[4]) +
    geom_text(aes(label = number), nudge_y = 50, size = 2) +
    scale_y_continuous(breaks = seq(0, 1300, 200), labels = seq(0, 1300, 200),
                       name = "frequency") +
    scale_x_discrete(name = "circRNA isoform number") +
    theme_classic() +
    mytheme
ggsave(paste0(figdir, "Figure2/isoform.pdf"), width = 75, height = 50, units = "mm")


########################################################################################################
# expression
gene_exp %>%
    pivot_longer(cols = WTM1:WTY3, 
                 names_to = "sample",
                 values_to = "CPG") %>%
    mutate(condition = ifelse(grepl("M", sample), "mycelium", "yeast")) %>%
    group_by(gene_id, condition) %>%
    summarise(mean = mean(CPG)/1000) %>%
    mutate(type = ifelse(grepl("TM", gene_id), "gene", "circRNA")) %>%
    ggplot(., aes(x = mean, fill = condition, color = condition)) +
        geom_vline(xintercept = c(0.01, 0.1, 1, 10, 100, 1000, 10000),
                   linetype = "dashed", color = "grey80", size = 0.15) +
        geom_density(size = 0.25) +
        scale_x_log10() +
        labs(x = "expression level (Counts Per Million)") +
        scale_fill_manual(values = alpha(col_nejm[1:2], 0.15)) +
        scale_color_manual(values = col_nejm[1:2]) +
        facet_grid(~type) +
        theme_bw() +
        mytheme +
        theme(
            legend.position = c(0.85,0.85),
            legend.background = element_blank(),
            legend.title = element_blank(),
            panel.grid = element_blank(),
            strip.background = element_rect(fill = NA)
        )
ggsave(paste0(figdir, "Figure2/CPM.pdf"), width = 95, height = 50, units = "mm")

# ratio
circ_ratio_long %>%
    group_by(circ_id, condition) %>%
    summarise(mean_ratio = mean(ratio)) %>%
    # pivot_wider(names_from = condition, values_from = mean_ratio) %>%
    # left_join(circ_type) %>%
    # ggplot(aes(x = Mycelium, y = Yeast, color = circ_type)) +
    ggplot(aes(x = mean_ratio, fill = condition, color = condition)) +
    geom_vline(xintercept = 0.1,
               linetype = "dashed", color = "grey80", size = 0.15) +
    # geom_point(size = 1.2, alpha = 0.1, color = "grey66") +
    # geom_point(size = 1, alpha = 0.5) +
    #labs(x = "ratio of mycelium", y = "ratio of yeast") +
    geom_density(size = 0.25) +
    scale_x_log10() +
    labs(x = "circular/linear ratio") +
    scale_fill_manual(values = alpha(col_nejm[1:2], 0.15)) +
    scale_color_manual(values = col_nejm[1:2]) +
    facet_grid(~condition) +
    theme_bw() +
    mytheme +
    theme(
        legend.position = c(0.85,0.85),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = NA)
    )
ggsave(paste0(figdir, "Figure2/ratio_2.pdf"), width = 95, height = 50, units = "mm")


########################################################################################################
# pm1_01:2332393|2332684 TM010856.00: EF-1 alpha subunit
# pm1_02:869646-869895 TM020307

# yeast pm1_07:211671|212431
# mycelium pm1_01:6114499|6115238

# cpm and ratio
circ_ratio_long_mean = circ_ratio_long %>%
    mutate(condition = ifelse(condition == "Mycelium", "mycelium", "yeast")) %>%
    group_by(circ_id, condition) %>%
    summarise(mean_ratio = mean(ratio))

circ_cpm_long_mean = circ_cpg_long %>%
    mutate(cpm = cpg / 1000) %>%
    group_by(circ_id, condition) %>%
    summarise(mean_cpm = mean(cpm))
    
circ_ratio_cpm_mean = left_join(circ_cpm_long_mean, circ_ratio_long_mean)

circ_ratio_cpm_mean %>%
    filter(condition == "mycelium") %>%
    arrange(-mean_ratio, -mean_cpm)

circ_ratio_cpm_mean = circ_ratio_cpm_mean %>%
    mutate(spe = ifelse(circ_id == "pm1_02:1328636|1328808", "circTM020485", ""),
           spe = ifelse(circ_id == "pm1_01:6114499|6115238", "circF1ATPase-a", spe),
           sz = ifelse(circ_id == "pm1_02:1328636|1328808", 0.75, 0.5),
           sz = ifelse(circ_id == "pm1_01:6114499|6115238", 0.75, sz))

circ_ratio_cpm_mean %>% 
    filter(condition == "yeast") %>%
    ggplot(aes(x = mean_ratio, y = mean_cpm, color = condition)) +
    geom_vline(xintercept = 0.1,
               linetype = "dashed", color = "grey60", size = 0.15) +
    geom_hline(yintercept = 0.1,
               linetype = "dashed", color = "grey60", size = 0.15) +
    geom_point(aes(size = sz)) +
    geom_text(aes(label = spe), color = "black", size= 3) +
    scale_x_log10(name = "circular/linear ratio") +
    scale_y_log10(name = "", sec.axis = sec_axis(~.*1 , name = "")) +
    labs(y = "expression level (Counts Per Million) under mycelium") +
    scale_color_manual(values = alpha(col_nejm[1], 0.5)) +
    theme_bw() +
    mytheme +
    theme(
        legend.position = "none",
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = NA),
    )
ggsave(paste0(figdir, "Figure2/mean_ratio_cpm_yeast.pdf"), width = 90, height = 70, units = "mm")

circ_ratio_cpm_mean %>% 
    filter(condition == "yeast") %>%
    ggplot(aes(x = mean_ratio, color = condition, fill = condition)) +
    geom_vline(xintercept = 0.1,
               linetype = "dashed", color = "grey60", size = 0.15) +
    geom_density(size = 0.25) +
    scale_x_log10() +
    labs(x = "", y = "") +
    scale_color_manual(values = col_nejm[1]) +
    scale_fill_manual(values = alpha(col_nejm[1], 0.25)) +
    theme_bw() +
    mytheme +
    theme(
        legend.position = "none",
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = NA),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
    )

ggsave(paste0(figdir, "Figure2/mean_ratio_density_yeast.pdf"), width = 75, height = 25, units = "mm")   
    

circ_ratio_cpm_mean %>% 
    filter(condition == "yeast") %>%
    ggplot(aes(x = mean_cpm, color = condition, fill = condition)) +
    geom_vline(xintercept = 0.1,
               linetype = "dashed", color = "grey60", size = 0.15) +
    geom_density(size = 0.25) +
    scale_x_log10() +
    labs(x = "", y = "") +
    scale_color_manual(values = col_nejm[1]) +
    scale_fill_manual(values = alpha(col_nejm[1], 0.25)) +
    theme_bw() +
    mytheme +
    theme(
        legend.position = "none",
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = NA),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
    )
ggsave(paste0(figdir, "Figure2/mean_cpm_density_yeast.pdf"), width = 71, height = 25, units = "mm")   



##
p = circ_ratio_cpm %>% 
    filter(condition == "mycelium") %>%
    ggplot(aes(x = ratio, y = cpm, color = condition)) +
    geom_vline(xintercept = 0.1,
               linetype = "dashed", color = "grey60", size = 0.15) +
    geom_hline(yintercept = 0.1,
               linetype = "dashed", color = "grey60", size = 0.15) +
    geom_point(size = 0.5) +
    scale_x_log10(name = "circular/linear ratio") +
    scale_y_log10() +
    #scale_y_log10(name = "", sec.axis = sec_axis(~.*1 , name = "")) +
    #labs(y = "expression level (Counts Per Million) under mycelium") +
    scale_color_manual(values = alpha(col_nejm[2], 0.5)) +
    theme_bw() +
    mytheme +
    theme(
        legend.position = "none",
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = NA),
    )

p2 = ggMarginal(p, type="density", 
                color = col_nejm[2], 
                fill = alpha(col_nejm[2], 0.5))

# Show only marginal plot for x axis
p3 <- ggMarginal(p2, margins = 'y', type="density", 
                 color = col_nejm[2], 
                 fill = alpha(col_nejm[2], 0.5))
