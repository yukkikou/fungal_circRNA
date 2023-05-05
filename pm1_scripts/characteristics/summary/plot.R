#!/usr/bin/Rscript

# signature
## Author: yuki
## Date: 2023.03.04

## environment
workdir = "D://_CircRNA/results/10_revise/20230201/"
rltdir="D://_CircRNA/results/"
figdir = "D://_CircRNA//figure//_panel//revise//20230210/"
set.seed(24)
setwd(workdir)

## packages
source("D://_Scripts/R/_LIBRARY.R")
library(paletteer)
library(rlang)
library(ggpubr)
library(patchwork)
library(ggthemes)

# colors and themes
source("D://_Scripts/R/_COLOR/colors.R")
source("D://_Scripts/R/_COLOR/ggsci_color_theme.R")

# functions
geneModCnt = function(geneMod, mod, con, type){
    mod = paste(con, str_split(mod, "ME", simplify = T)[,2], sep = "_")
    if (type == "gene"){
        geneMod %>% filter(module == mod, grepl("TM",feat_id)) %>%
            count()
    } else if (type == "circ"){
        geneMod %>% filter(module == mod, grepl("pm", feat_id)) %>%
            count()
    } else {
        stop("Error")
    }
}

cytoExt10 = function(circEdges, gene){
    circEdges %>% filter((fromNode == gene)|(toNode == gene)) %>%
        head(10)
}


## Result 1
# Figure 1B and FigureS1: F1 score in mouse
mouse_ciri_id = read_tsv(paste0(rltdir, "/11_mouse/3_circAltas/CIRI2.id"), col_names = "mus_circ_id")
mouse_long_id = read_tsv(paste0(rltdir, "/11_mouse/3_circAltas/CIRIlong.id"), col_names = "mus_circ_id")
mouse_inter_id = read_tsv(paste0(rltdir, "/11_mouse/3_circAltas/CIRI2_CIRIlong.id"), col_names = "mus_circ_id")
mouse_merge_id = read_tsv(paste0(rltdir, "/11_mouse/1_NanoIllu/merge_all.id"), col_names = "mus_circ_id")
mouse_alt_id = read_tsv(paste0(rltdir, "/11_mouse/3_circAltas/mouse.clean.bed"), 
                        col_names = c("chr", "ss", "es", "mus_circ_id", "score", "strand"))

mouse_alt_number = dim(mouse_alt_id)[1]

overlap_method = c("Illumina","Nanopore", "Intersection", "Recorrection")
identication_data = data.frame(
    "method" = overlap_method,
    "find_number" = c(dim(mouse_ciri_id)[1], dim(mouse_long_id)[1], dim(mouse_inter_id)[1], dim(mouse_merge_id)[1]),
    "overlap_number" = c(table(mouse_ciri_id$mus_circ_id %in% mouse_alt_id$mus_circ_id)[2],
                         table(mouse_long_id$mus_circ_id %in% mouse_alt_id$mus_circ_id)[2],
                         table(mouse_inter_id$mus_circ_id %in% mouse_alt_id$mus_circ_id)[2],
                         table(mouse_merge_id$mus_circ_id %in% mouse_alt_id$mus_circ_id)[2])
) %>% mutate(
    "TP" = overlap_number/find_number,
    "FP" = (find_number - overlap_number)/find_number,
    "FN" = (mouse_alt_number - overlap_number)/mouse_alt_number,
    "F1score" = round(2*TP/(2*TP + FP + FN), 2),
)

identication_data_long = identication_data %>% pivot_longer(cols = find_number:F1score,
                                   names_to = "args",
                                   values_to = "value") %>%
    mutate(args = factor(args, levels = c("TP", "FP", "FN", "F1score", "overlap_number", "find_number")),
           method = factor(method, levels = overlap_method),
           args_type = ifelse(grepl("number", args), "Number", "statistics")) 

identication_data_long %>% 
    filter(args_type == "Number") %>%
    ggplot(aes(x = method, y = value, fill = method)) +
    facet_grid(cols = vars(args),
               scales = "free") +
        geom_col(position = position_dodge(), width = 0.18) +
        geom_text(aes(label = round(value, 2)), size = 2,
                  color = "grey5", nudge_y = 1000) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 18000)) +
        scale_x_discrete(labels = paste0("M",1:4)) +
        scale_fill_manual(values = col_nejm) +
        labs(x = "", y = "") +
        theme_bw() +
        mytheme +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            strip.background = element_rect(color="grey33", fill = "grey95"),
        )
ggsave(paste0(figdir, "Figure1/recor_number.pdf"), width = 90, height = 50, units = "mm")

identication_data_long %>% 
    filter(args == "F1score") %>%
    ggplot(aes(x = method, y = value, fill = method)) +
    facet_grid(cols = vars(args),
               scales = "free") +
    geom_col(position = position_dodge(), width = 0.18) +
    geom_text(aes(label = round(value, 2)), size = 2,
              color = "grey5", nudge_y = 0.05) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.8)) +
    scale_x_discrete(labels = paste0("M",1:4)) +
    scale_fill_manual(values = col_nejm) +
    labs(x = "", y = "") +
    theme_bw() +
    mytheme +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_rect(color="grey33", fill = "grey95"),
    )
ggsave(paste0(figdir, "Figure1/recor_f1score.pdf"), width = 45, height = 50, units = "mm")

identication_data_long %>% 
    filter(args != "F1score", args_type == "statistics") %>%
    ggplot(aes(x = method, y = value, fill = method)) +
    facet_grid(cols = vars(args),
               scales = "free") +
    geom_col(position = position_dodge(), width = 0.25) +
    geom_text(aes(label = round(value, 2)), size = 2,
              color = "grey5", nudge_y = 0.05) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.2)) +
    scale_x_discrete(labels = paste0("M",1:4)) +
    scale_fill_manual(values = col_nejm) +
    labs(x = "", y = "") +
    theme_bw() +
    mytheme +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_rect(color="grey33", fill = "grey95"),
        panel.spacing = unit(1.5, "lines")
    )
ggsave(paste0(figdir, "FigS1/recor_score.pdf"), width = 150, height = 50, units = "mm")

# Figure 1C
mycelium_circ_id = read_tsv(paste0(rltdir, "8_circ_WT_merge/1_number/mycelium_all.vertical.circid"), col_names = "circ_id")
yeast_circ_id = read_tsv(paste0(rltdir, "8_circ_WT_merge/1_number/yeast_all.vertical.circid"), col_names = "circ_id")
share_circ_id = inner_join(mycelium_circ_id, yeast_circ_id)

circ_type = circ_info %>%
    mutate(circ_type = ifelse(circ_id_vertical %in% share_circ_id$circ_id, "share", 
                              ifelse(circ_id_vertical %in% mycelium_circ_id$circ_id, "mycelium", "yeast")),
           nano_supported = ifelse(is.na(len_same), "no", "yes")) %>%
    select(circ_id_vertical, nano_supported, circ_type) %>%
    rename(circ_id = circ_id_vertical)

circ_type[,-1] %>% table() %>% as.data.frame() %>%
    ggplot(aes(x = circ_type, y = Freq, fill = nano_supported)) +
        geom_hline(yintercept = c(500, 1500, 2500), 
                   linetype = "dashed", color = "grey60",
                   size = 0.35) +
        geom_col(position = position_stack(), width = 0.18) +
        scale_y_continuous(breaks = seq(0,3000,500), name = "") +
        scale_x_discrete(name = "", labels = c("M-sp", "Share", "Y-sp")) +
        scale_fill_manual(values = col_nejm[c(2,1)]) +
        theme_bw() +
        mytheme +
        theme(
            legend.position = c(0.32,0.73),
            panel.grid = element_blank(),
            legend.title = element_text(size = 7),
            legend.background = element_rect(fill = "NA")
        )
ggsave(paste0(figdir, "Figure1/nanosupp_number.pdf"), width = 50, height = 50, units = "mm")

# circ_cpm_mean = circ_cpg_long %>% 
#     group_by(circ_id, condition) %>%
#     summarise(log10CPG = log10(mean(cpg)+1)) %>%
#     left_join(circ_type)
# 
# circ_cpm_mean %>% 
#     ggplot(aes(x = circ_type, y = log10CPG, color = nano_supported)) +
#     geom_boxplot(width = 0.4) +
#     facet_grid(rows = vars(condition)) +
#     scale_color_manual(values = alpha(col_nejm[c(2, 1)], 0.5)) +
#     labs(x = "Types of circRNAs", y = "log10(CPG+1)") +
#     theme_bw() +
#     mytheme +
#     theme(
#         legend.position = c(0.3, 0.9),
#         strip.background = element_rect(color="grey33", fill = "grey95", size = 0)
#     )
# ggsave(paste0(figdir, "Figure1/nanosupp.pdf"), width = 65, height = 100, units = "mm")

## Result2
# Figure2A
gene_length = read_tsv(paste0(rltdir, "8_circ_WT_merge/2_length/pm1.gene.length"), 
                       col_names = c("chr", "len", "gene_id")) %>%
    mutate("feat_id" = gene_id,
           "feat_type" = "gene") %>%
    select("feat_id", "len", "gene_id", "feat_type")
circ_length = circ_info %>% select("circ_id_vertical", "len", "gene_id") %>%
    mutate("feat_type" = "circRNA") %>%
    filter(!is.na(len)) %>%
    rename(feat_id = circ_id_vertical) %>%
    select("feat_id", "len", "gene_id", "feat_type")

circ_gene_len = inner_join(circ_length, gene_length,
          by = c("gene_id" = "gene_id"),
          suffix = c(".circ", ".gene")) %>%
    mutate(log10circ = log10(len.circ),
           log10gene = log10(len.gene))


rbind(gene_length, circ_length) %>%
    ggplot(aes(x = len, color = feat_type, fill = feat_type)) +
        geom_density(size = 0.25) +
        scale_color_manual(values = col_nejm[3:4]) +
        scale_fill_manual(values = alpha(col_nejm[3:4], 0.3)) +
        scale_x_log10() +
        labs(x  = "length (nt)") +
        theme_bw() +
        mytheme +
        theme(
            legend.position = c(0.8,0.85),
            legend.title = element_blank(),
            legend.background = element_blank()
        )
ggsave(paste0(figdir, "Figure2/circ_gene_length1.pdf"), width = 55, height = 48, units = "mm")
ggplot(circ_gene_len, aes(x = log10gene, y = log10circ)) +
    geom_smooth(method = "lm", size = 0.35, color = col_nejm[1]) +
    labs(x = "log10(length of gene)", y = "log10(length of circ)") +
    theme_bw() +
    mytheme
ggsave(paste0(figdir, "Figure2/circ_gene_length2.pdf"), width = 50, height = 48, units = "mm")    
cor.test(circ_gene_len$log10circ, circ_gene_len$log10gene,method = "spearman")

# Figure2B
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
    geom_boxplot(width  = 0.5, outlier.alpha = 0.5, outlier.size = 0.5) +
    scale_color_manual(values = rev(col_nejm)[c(4,3,6)]) + 
    labs(x = "circ type") +
    theme_bw() +
    mytheme +
    theme(
        legend.position = c(0.8,0.8),
        legend.background = element_rect(fill = NA),
        legend.title = element_blank()
    )
ggsave(paste0(figdir, "Figure2/circ_length.pdf"), width = 55, height = 48, units = "mm")

t.test(filter(circ_length_compare, type != "NSD")$length,
       filter(circ_length_compare, type == "NSD")$length)

# Figure2C
# linear isoform and circular isoform
gene_transnumber = read_tsv("D://_IGV/Talaromyces_marneffei/Annotation/gene_linear.number", col_names = c("linear_number", "gene_id"))
gene_circnumber = table(circ_info$gene_id) %>% as.data.frame() %>% `colnames<-`(c("gene_id", "circular_number"))
gene_isonumber = inner_join(gene_transnumber, gene_circnumber) %>% select(gene_id, linear_number, circular_number) %>%
    mutate(all_number = linear_number + circular_number,
           circular_per = round(circular_number/all_number*100,1))

cor.test(gene_isonumber$circular_number, gene_isonumber$linear_number, method = "spearman")

ggplot(gene_isonumber, aes(x = circular_number, y = linear_number)) +
    geom_point(size = 0.5, color = alpha(col_nejm[4], 0.3)) +
    geom_smooth(color = col_nejm[2], size = 0.5) +
    # stat_regline_equation(label.x = 22, label.y = 90, aes(label = ..eq.label..)) +
    # stat_regline_equation(label.x = 22, label.y = 80, aes(label = ..rr.label..)) +
    labs(x = "circular isoform number", y = "linear isoform number") +
    theme_bw() +
    mytheme
ggsave(paste0(figdir, "Figure2/circ_linear.pdf"), width = 70, height = 45, units = "mm")

# Figure2D
circ_ratio_long %>%
    group_by(circ_id, condition) %>%
    summarise(mean_ratio = mean(ratio)) %>%
    pivot_wider(names_from = condition, values_from = mean_ratio) %>%
    left_join(circ_type) %>%
    ggplot(aes(x = Mycelium, y = Yeast, color = circ_type)) +
        geom_point(size = 2, alpha = 0.1, color = "grey66") +
        geom_point(size = 1.8, alpha = 0.5) +
        labs(x = "ratio of mycelium", y = "ratio of yeast") +
        scale_color_manual(values = col_nejm) +
        theme_bw() +
        mytheme +
        theme(
            legend.title = element_blank(),
            legend.position = "top"
        )
ggsave(paste0(figdir, "Figure2/circ_rario.pdf"), width = 70, height = 75, units = "mm")

## Result3
# Figure3B
marker_number_map = data.frame(
    circ_marker = c("intergenic", "single_exon_over", "single_exon_match", "single_exon_internal", 
                    "multi_single_exons", 
                    "multi_single_internal_exon_exact_match", "multi_single_edge_exon_exact_match", 
                    "multi_single_edge_over",
                    "multi_intron_start", "multi_intron_end", "multi_intron_start_end",
                    "multi_single_intron_match", "multi_single_intron_internal", 
                    "multi_exon_exon", "multi_exon_intron_exons",
                    "ambiguous", "antisense", "read_through"),
    circ_code = c(1:6, 6, 7:12, 14, 13, 15:17)
)

gene_exon_type = rbind(read_tsv(paste0(rltdir, "8_circ_WT_merge/5_structure/single_exon.gene"), col_names = "gene_id") %>%
          mutate(gene_type = "single_exon"), 
    read_tsv(paste0(rltdir, "8_circ_WT_merge/5_structure/multi_exon.gene"), col_names = "gene_id") %>%
        mutate(gene_type = "multi_exon"))

circ_marker_table = read_tsv(paste0(rltdir, "8_circ_WT_merge/5_structure/circ_marker_table.tsv")) %>%
    mutate(gene_id = paste0(str_split(transcript_id, pattern = "\\.0[1-9]", simplify = T)[,1], ".00")) %>%
    left_join(gene_exon_type)

duplicated_circ = circ_marker_table[duplicated(circ_marker_table$circ_id),]
circ_marker_table %>% filter(circ_id %in% duplicated_circ$circ_id) %>% arrange(circ_id) # read though

circ_marker_table_clean = circ_marker_table %>% filter(!(circ_id %in% duplicated_circ$circ_id))
table(circ_marker_table_clean$gene_type)
table(circ_marker_table_clean[,c(4,6)])

circ_marker_table_info = circ_info %>% filter(!is.na(len)) %>% select(circ_id_hyphen, gene_id, circ_type) %>%
    left_join(circ_marker_table_clean, by = c("circ_id_hyphen" = "circ_id", "gene_id" = "gene_id")) %>%
    mutate(circ_type = ifelse(circ_id_hyphen %in% duplicated_circ, "read_though", circ_type),
           gene_id = ifelse(is.na(gene_id), circ_id_hyphen, gene_id),
           gene_type = ifelse(grepl("pm1", gene_id), "intergenic", gene_type),
           circ_type = ifelse(grepl("pm1", gene_id), "intergenic", circ_type),
           circ_marker = ifelse(grepl("pm1", gene_id), "intergenic", circ_marker),
           circ_marker = ifelse(circ_type == "antisense", "antisense", circ_marker),
           circ_marker = ifelse(is.na(circ_marker), "ambiguous", circ_marker)) %>%
    left_join(marker_number_map)

table(circ_marker_table_info$circ_code) %>%
    as.data.frame() %>%
    `colnames<-`(c("circ_code", "number")) %>%
    mutate(circ_code = factor(circ_code, levels = 1:16)) %>%
    ggplot(aes(x = circ_code, y = number, fill = as.character(c(1, rep(2, 3), rep(3, 8), rep(4, 2))))) +
        geom_hline(yintercept = c(200,400,600), color = "grey66", size = 0.25, linetype = "dashed") +
        geom_col(width = 0.3) +
        geom_text(aes(label = number), nudge_y = 40, color = "grey30", size = 2.5) +
        scale_fill_manual(values = col_nejm) +
        scale_y_continuous(expand = c(0,0), limits = c(0, 750), name = "") +
        scale_x_discrete(name = "circRNA type") +
        theme_classic() +
        #theme_bw() +
        mytheme +
        theme(
            legend.position = "none"
        )

ggsave(paste0(figdir, "Figure3/number_type_clean.pdf"), width = 100, height = 50, units = "mm")

# Figure3C
intron_number_length = read_tsv(paste0(rltdir, "8_circ_WT_merge/5_structure/intron_number_t01.bed"), 
                                col_names = c("chr", "ss", "es", "intron_id", "intron_length", "strand"))
circ_intron_feat = read_tsv(paste0(rltdir, "8_circ_WT_merge/5_structure/circ_feat_intron.tsv"), col_names = F)[,c(3,8)] %>%
    `colnames<-`(c("circ_id", "intron_id"))

gene_circ_intron = full_join(circ_intron_feat, intron_number_length) %>%
    mutate(feat_type = ifelse(is.na(circ_id), "gene", "circRNA"))

ggplot(gene_circ_intron, aes(x = log10(intron_length), color = feat_type)) +
    geom_density(size = 0.5) +
    scale_color_manual(values = col_nejm) +
    scale_y_continuous(name = "") +
    scale_x_continuous(name = "log10(length of intron)") +
    theme_classic() +
    mytheme +
    theme(
        legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA,),
        panel.grid = element_blank()
    )

ggsave(paste0(figdir, "Figure3/intron_length1_clean.pdf"), width = 50, height = 40, units = "mm")

ggplot(gene_circ_intron, aes(x = feat_type, y = log10(intron_length), color = feat_type)) +
    geom_boxplot(width = 0.5, outlier.size = 0.5, outlier.alpha = 0.5) +
    scale_color_manual(values = col_nejm) +
    scale_y_continuous(name = "") +
    scale_x_discrete(name = "log10(length of intron)") +
    theme_classic() +
    mytheme +
    theme(
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA)
    )

ggsave(paste0(figdir, "Figure3/intron_length2_clean.pdf"), width = 50, height = 40, units = "mm")

wilcox.test(filter(gene_circ_intron, feat_type == "circRNA")$intron_length,
            filter(gene_circ_intron, feat_type == "gene")$intron_length)

# Figure 3D
circ_ssite = read_tsv(paste0(rltdir, "8_circ_WT_merge/4_spliceSite/all.motifRev.table"), 
                      col_names = c("motif", "count", "freq")) %>%
    mutate(type = "circRNA")
gene_ssite = read_tsv(paste0(rltdir, "8_circ_WT_merge/4_spliceSite/pm1_intron.motifRev.table"),
                      col_names = c("motif", "count", "freq")) %>%
    mutate(type = "gene",
           freq = -freq,
           count = -count)
rbind(circ_ssite, gene_ssite) %>%
    filter(motif != "AGGT") %>%
    filter(abs(freq) > 0.0005) %>%
    ggplot(aes(x = motif, y = freq*100, fill = type)) +
        geom_col(width = 0.5) +
        scale_x_discrete(name = "") +
        scale_y_continuous(name = "") +
        scale_fill_manual(values = col_jama[c(5,1)]) +
        theme_bw() +
        mytheme +
        theme(
            axis.text.x = element_text(angle = 45, size = 6),
            legend.position = c(0.8,0.8),
            legend.background = element_rect(fill = NA),
            panel.grid = element_blank()
        )
ggsave(paste0(figdir, "Figure3/motif_non.pdf"), width = 90, height = 60, units = "mm")

rbind(circ_ssite, gene_ssite) %>%
    filter(motif == "AGGT") %>%
    ggplot(aes(x = motif, y = freq*100, fill = type)) +
    geom_col(width = 0.5) +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Percentage") +
    scale_fill_manual(values = col_jama[c(5,1)]) +
    theme_bw() +
    mytheme +
    theme(
        axis.text.x = element_text(angle = 45, size = 6),
        legend.position = "none",
        legend.background = element_rect(fill = NA),
        panel.grid = element_blank()
    )
ggsave(paste0(figdir, "Figure3/motif_can.pdf"), width = 22, height = 60, units = "mm")

# Figure 3d
keep_motif = rbind(circ_ssite, gene_ssite) %>%
    filter(type == "circRNA") %>%
    filter(abs(freq) > 0.0005)

# tmp classification
circ_ssite_table = read_tsv(paste0(rltdir, "8_circ_WT_merge/4_spliceSite/all.motifRev.tsv"))
circ_ss_marker_clean = inner_join(circ_ssite_table, circ_marker_table_clean,
                                  by = c("circID" = "circ_id")) %>% 
    left_join(marker_number_map) %>% 
    mutate(motif_type = ifelse(motifRev == "AGGT", "canonical", "non-canonical"),
           border_type = ifelse(grepl("match", circ_marker), "match", "dismatch"),
           border_type = ifelse(circ_code %in% 13:14, "dismatch", border_type))

circ_ss_marker_clean = circ_ss_marker_clean %>% filter(!(circ_code %in% 13:14))

circ_ss_marker_clean %>%
    select(motif_type, border_type) %>%
    table() %>% chisq.test()

circ_ss_marker_clean %>%
    filter(motif_type == "non-canonical", border_type == "dismatch",
           motifRev == "AGAT") %>%
    select(circID, compose:border_type)

circ_ss_marker_clean %>%
    filter(motif_type == "non-canonical", border_type == "dismatch",
           motifRev == "AGGC") %>%
    select(circID, compose:border_type)

circ_ss_marker_clean_non = circ_ss_marker_clean %>%
    filter(motif_type == "non-canonical") %>%
    select(circID, motifRev, compose:border_type)

table(circ_ss_marker_clean_non[,c(4,9)])
table(circ_ss_marker_clean_non[,c(2,9)]) %>% 
    as.data.frame() %>%
    filter(motifRev %in% circ_ss_marker_clean$motifRev[table(circ_ss_marker_clean$motifRev)>4]) %>%
    ggplot(aes(x = motifRev, y = Freq, fill = border_type)) +
        geom_hline(yintercept = c(10, 20), linetype = "dashed", size = 0.5, color = "grey66") +
        geom_col(width = 0.5) +
        scale_fill_manual(values = col_nejm) +
        scale_x_discrete(name = "") +
        scale_y_continuous(name = "frequency", breaks = seq(0, 40, 5), labels = seq(0, 40, 5), expand = c(0, 0)) +
        coord_flip() +
        theme_classic() +
        mytheme +
            theme(
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = c(0.85, 0.85),
                #legend.background = element_blank()
            )

#ggsave(paste0(figdir, "Figure3/motif_non_ss.pdf"), width = 60, height = 100, units = "mm")
ggsave(paste0(figdir, "Figure3/motif_non_ss_filter.pdf"), width = 60, height = 60, units = "mm")

# Result4
motif_anno = data.frame(
    motif_struc = unique(c(hnRNP_motif, SF_motif)),
    motif_anno = c("hnRNP,SF","hnRNP","hnRNP,SF","hnRNP","hnRNP", rep("SF", 3))
) %>% left_join(motif_bind[grepl("HNRNP", motif_bind$gene_id),], by = c("motif_struc" = "motif")) %>%
    left_join(motif_bind[grepl("SF", motif_bind$gene_id),], by = c("motif_struc" = "motif"),
              suffix = c(".hnRNP", ".SF"))
human_pm1_map = read_tsv(paste0(rltdir, "8_circ_WT_merge/6_flanking/3_motif/3_homoGene/human_pm1.gene"),
                         col_names = c("human.gene", "human.alias", "pm1", "func", "homo")) %>%
    pivot_longer(cols = c(human.gene, human.alias),
                 names_to = "human",
                 values_to = "human_gene") %>%
    select(-human)
 
circ_flk_table = rbind(mutate(circ_flk_down, type = "down_stream"), mutate(circ_flk_up, type = "up_stream"))
circ_flk_table %>%
    filter(motif_struc %in% c(hnRNP_motif, SF_motif)) %>%
    left_join(motif_anno) %>%
    pivot_longer(cols = c(gene_id.hnRNP, gene_id.SF),
                 names_to = "human_gene_type",
                 values_to = "human_gene") %>%
    left_join(human_pm1_map) %>%
    filter(!is.na(human_gene)) %>%
    select(-human_gene_type) %>%
    write_tsv("D://_CircRNA/supplementary/20230210/SupplementaryTable1.tsv")

circ_flk_table_frq = circ_flk_table %>% group_by(type, motif_id, p) %>%
    summarise(count = n()) %>%
    mutate(precent = count/3891*100,
           axis = ifelse(type == "down_stream", precent, -precent),
           id = paste0(type, motif_id))

circ_flk_table_frq %>%
    filter(precent > 4.5) %>%
    ggplot(aes(x = id, y = axis, alpha = p)) +
        geom_col(width = 0.3) +
        scale_x_discrete(labels = NULL, name = "") +
        scale_y_continuous(name = "", labels = paste0(c(rev(seq(0,20,10)),seq(5,20,10)), "%")) +
        coord_flip()+
        theme_minimal() +
        mytheme +
        theme(
            #legend.position = "none",
            #legend.title = element_blank(),
            legend.background = element_rect(fill = NA),
            panel.grid = element_blank()
        )

ggsave(paste0(figdir, "Figure3/circ_flk_table_frq.pdf"), width = 90, height = 70, units = "mm")

## Result5
# figure4A
de_circ_id = circ_de %>% filter(DE != 0)

circ_bsj %>% filter(circ_id %in% de_circ_id$circ_id_vertical) %>%
    column_to_rownames("circ_id") %>%
    pheatmap::pheatmap(show_rownames = F, scale = "row",
                       colorRampPalette(c(col_nejm[2], "white",  col_nejm[1]))(50), 
                       treeheight_row = 15, treeheight_col = 15)

de_gene_id = read_tsv(paste0(rltdir, "8_circ_WT_merge/7_expression/1_linear/gene/gene_res_sig.tsv"))
gene_tpm = read_csv(paste0(rltdir, "8_circ_WT_merge/7_expression/2_circ/WT_merge_gene_TPM_matrix.csv")) %>%
    filter(gene_id %in% de_gene_id$gene_id)

gene_tpm %>% column_to_rownames("gene_id") %>%
    pheatmap::pheatmap(show_rownames = F, scale = "row",
                       colorRampPalette(c(col_nejm[2], "white",  col_nejm[1]))(50), 
                       treeheight_row = 15, treeheight_col = 15)


# figure4B
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
    scale_y_discrete(labels = \(x) {str_wrap(x, width = 25)}, name = "") +
    scale_color_continuous(type = "viridis") +
    theme_bw() +
    mytheme +
    theme(
        axis.text.y = element_text(size = 7),
    )
ggsave(paste0(figdir, "Figure4/go_circ.pdf"), width = 80, height = 55, units = "mm")

de_gene_go %>% slice(1:5) %>%
    ggplot(., aes(x = FE,y = Description, color = p.adjust)) +
    geom_point() +
    scale_x_continuous(limits = c(0, 0.12), name = "fold enrichment") +
    scale_y_discrete(labels = \(x) {str_wrap(x, width = 25)}, name = "") +
    scale_color_continuous(type = "viridis") +
    theme_bw() +
    mytheme +
    theme(
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(angle = 45)
    )
ggsave(paste0(figdir, "Figure4/go_gene.pdf"), width = 85, height = 63, units = "mm")

# figure4C
de_gene = read_tsv(paste0(rltdir, "8_circ_WT_merge/7_expression/1_linear/gene/gene_res.tsv"))


# volcano plot
de_gene_res_anno = de_gene %>% inner_join(circ_info, by = c("gene_id" = "gene_id")) %>%
    select(gene_id, log2FoldChange, padj, circ_id_vertical) %>%
    dplyr::rename(geneFC = log2FoldChange,
                  geneFDR = padj)
de_circ_res_anno = circ_de %>% 
    select(circ_id_vertical, logFC, FDR) %>%
    dplyr::rename(circFC = logFC,
                  circFDR = FDR)
de_gene_circ_table = left_join(de_gene_res_anno, de_circ_res_anno) %>%
    mutate(type = NA)

for (i in 1:dim(de_gene_circ_table)[1]){
    tmp_geneFC = abs(de_gene_circ_table[i, "geneFC"])>1
    tmp_circFC = abs(de_gene_circ_table[i, "circFC"])>1
    if(tmp_geneFC & tmp_circFC){
        de_gene_circ_table[i, "type"] = "both"
    } else if (tmp_geneFC) {
        de_gene_circ_table[i, "type"] = "gene"
    } else if (tmp_circFC) {
        de_gene_circ_table[i, "type"] = "circ"
    } else {
        de_gene_circ_table[i, "type"] = "none"
    }
}

de_gene_circ_table %>%
    ggplot(aes(x = geneFC, y = circFC, color = type)) +
    geom_hline(yintercept = 0, size = 0.5, color = "snow2") +
    geom_vline(xintercept = 0, size = 0.5, color = "snow2") +
    geom_abline(slope = 1, size = 0.5, color = "snow2", linetype = "dashed") +
    #geom_point(size = 2.2, alpha = 0.2, color = "grey60") +
    geom_point(size = 1.1, alpha = 0.65) +
    xlim(c(-11,11)) +
    ylim(c(-11,11)) +
    scale_color_manual(values = c(col8[c(5,6,7)], "grey60")) +
    labs(x = "log2(FoldChange) of gene", y = "log2(FoldChange) of circRNA") +
    theme_bw() +
    mytheme +
    theme(
        legend.position = "bottom",
        panel.grid=element_blank(),
        legend.text = element_text(size = 8)
    )
ggsave(paste0(figdir, "Figure4/gene_circ_valcano.pdf"), width = 63, height = 70, units = "mm")

# figure4D
circ_clasf_map %>% select(circ_id_vertical, log2FoldChange.gene, logFC, log2FoldChange.ratio, classification) %>%
    rename(gene = log2FoldChange.gene, 
           circRNA = logFC, 
           ratio = log2FoldChange.ratio) %>%
    pivot_longer(cols = gene:ratio,
                 names_to = "DE_type",
                 values_to = "log2FC") %>%
    mutate(DE_type = factor(DE_type, levels = c("circRNA", "ratio", "gene"))) %>%
    ggplot(aes(x = DE_type, y = circ_id_vertical, fill = log2FC)) +
        geom_tile() +
        scale_fill_gradient2(low = "#4E79A7", mid = "white", high = "#B07AA1") +
        scale_x_discrete(name = "", labels = c("Circ","Ratio","Gene")) +
        scale_y_discrete(name = "gene/circRNA ID", label = NULL) +
        coord_flip() +
        facet_grid(~classification, space = "free", scales = "free") +
        theme_minimal() +
        mytheme +
        theme(
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            strip.background.x = element_rect(
                color = NA, fill="grey80", size=0.5, linetype="solid"
            ),
            axis.title.x = element_text(size = 7),
            axis.text.y = element_text(size = 6)

        )
ggsave(paste0(figdir, "Figure4/gene_circ_ratio.pdf"), width = 140, height = 55, units = "mm")

## Result5
# qPCR PLOTS are in qPCR.R


## Result6 
# Figure 6A (heatmap)
# Figure 6B (circular plot)
 
# GO enrichment
grep("w_de", ls(), value = T)
source("D://_CircRNA/scripts/Enrichment/s00_goFunction.R")    

all_disp = data.frame(matrix(ncol = 5, nrow = 0))
for (i in c("m", "y")){
    for (j in c("k", "o")){
        tmp = deExact(get(paste0(i, j,"w_de")), a = "all")
        tmp = go_result(tmp) %>% filter(pvalue < 0.01) %>%
            mutate(group = paste0(i, j, "w"))
        assign(paste0(i, j ,"w_go"), tmp)
        print(tmp)
        print(dim(tmp))
        all_disp = rbind(all_disp, tmp[1:10,])
    }
}

all_disp_direc = data.frame(matrix(ncol = 6, nrow = 0))
for (i in c("m", "y")){
    for (j in c("k", "o")){
        tmp = deExact(get(paste0(i, j,"w_de")), a = "all")
        tmp = go_result(tmp) %>% filter(pvalue < 0.01) %>%
            mutate(group = paste0(i, j, "w"))
        assign(paste0(i, j ,"w_go"), tmp)
        print(tmp)
        print(dim(tmp))
        all_disp = rbind(all_disp, tmp[1:10,])
    }
}

#all_disp %>% write_tsv("./../../../supplementary/20230210/SupplementaryTable2.tsv")
#all_disp %>% write_csv("../../../supplementary/20230210/SupplementaryTable2.csv")
#all_disp = read_csv("../../../supplementary/20230210/SupplementaryTable2.csv")

all_de_disp = all_disp[(duplicated(all_disp$ID)),]
all_disp %>%
    mutate(share = ifelse(ID %in% all_de_disp$ID, "yes", "no"),
           EF = DOSE::parse_ratio(GeneRatio)) %>%
    arrange(desc(share)) %>%
    mutate(Description = factor(Description, levels = unique(Description)),
           ID = factor(ID, levels = unique(ID))) %>%
 #   filter(share == "yes") %>%
    ggplot(aes(x = pvalue, y = Description)) +
    geom_point(aes(size = EF, color = group, shape = share)) +
    scale_x_continuous(name = "P value") +
    scale_y_discrete(name = "", labels = \(x) {str_wrap(x, width = 45)}) +
    scale_color_manual(values = col_nejm) +
    theme_minimal() +
    mytheme +
    theme(legend.position = "right",
          legend.background = element_blank())

ggsave(paste0(figdir, "Figure6/mutant_go_label_share.pdf"), width = 170, height = 180, units = "mm")
ggsave(paste0(figdir, "FigS2/mutant_go_label_share.pdf"), width = 170, height = 190, units = "mm")

 # WGCNA module

all_geneMod = rbind(read_tsv(paste0(rltdir, "9_WGCNA/2_WGCNA/mycelium/mycelium_geneMod.tsv")),
                    read_tsv(paste0(rltdir, "9_WGCNA/2_WGCNA/yeast/yeast_geneMod.tsv")))

myce_net_cor = read_tsv(paste0(rltdir, "9_WGCNA/2_WGCNA/mycelium/mycelium_module_train_cor_pvalue.tsv")) %>%
    filter(abs(mutant.cor) > 0.5) %>% 
    mutate(type = "mycelium",
           module_id = paste0("M-Mod", toupper(substr(str_split(module, "ME", simplify = T)[,2], 1, 1))))

for (i in 1:nrow(myce_net_cor)){
    myce_net_cor[i, "gene_mod_num"] = geneModCnt(all_geneMod, myce_net_cor[i, "module"], "mycelium", "gene")$n
    myce_net_cor[i, "circ_mod_num"] = geneModCnt(all_geneMod, myce_net_cor[i, "module"], "mycelium", "circ")$n
}

yst_net_cor = read_tsv(paste0(rltdir, "9_WGCNA/2_WGCNA/yeast/yeast_module_train_cor_pvalue.tsv")) %>%
    filter(abs(mutant.cor) > 0.5) %>%
    mutate(type = "yeast",
           module_id = paste0("Y-Mod", toupper(substr(str_split(module, "ME", simplify = T)[,2], 1, 1))))

for (i in 1:nrow(yst_net_cor)){
    yst_net_cor[i, "gene_mod_num"] = geneModCnt(all_geneMod, yst_net_cor[i, "module"], "yeast", "gene")$n    
    yst_net_cor[i, "circ_mod_num"] = geneModCnt(all_geneMod, yst_net_cor[i, "module"], "yeast", "circ")$n    
}

all_net_cor = rbind(myce_net_cor, yst_net_cor)[,c(2, 13:15)] %>%
    pivot_longer(cols = gene_mod_num:circ_mod_num, 
                 names_to = "num_type", values_to = "num") %>%
    mutate(num_y = ifelse(num_type == "circ_mod_num", -num, num))

ggplot(all_net_cor, aes(x = module_id, y = num_y, color = mutant.cor)) +
    geom_hline(yintercept = seq(0, 800, 300), size = 0.35, linetype = "dashed", color = "grey80") +
    geom_segment(data = filter(all_net_cor, num_type == "gene_mod_num"),
                 mapping = aes(x=module_id, xend=module_id, y=0, yend=num_y), size = 0.5) +
    geom_segment(data = filter(all_net_cor, num_type == "circ_mod_num"),
                 mapping = aes(x=module_id, xend=module_id, yend=0, y=num_y), size = 0.5) +
    geom_text(data = filter(all_net_cor, num_type == "gene_mod_num"),
              mapping = aes(label = num), nudge_y = 60, size = 2.5) +
    geom_text(data = filter(all_net_cor, num_type == "circ_mod_num"),
              mapping = aes(label = num), nudge_y = -50, size = 2.5) +
    geom_point(size=1.5) +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "circRNA/gene number", breaks = seq(-100, 1000, 100),
                       labels = c(100, seq(0, 1000, 100))) +
    scale_color_gradient2(name = "correlation") +
    coord_flip() +
    theme_classic() +
    mytheme +
    theme(
        legend.position = c(0.9,0.4),
        axis.text.x = element_text(angle = 45),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()
    )

ggsave(paste0(figdir, "Figure6/net_module.pdf"), width = 85, height = 60, units = "mm")

# WGCNA enrichment
wcgnsdir= paste0(rltdir, "9_WGCNA/2_WGCNA")
flist = list.files(path = wcgnsdir, pattern = "*_enrich.tsv", recursive = T, full.names = T)
all_mod_go = data.frame(matrix(ncol =10, nrow = 0))

for (i in flist){
    var = substr(str_split(str_split(i, "/", simplify = T)[,8], "_", simplify = T)[,2], 1, 1) %>% toupper()
    con = substr(str_split(str_split(i, "/", simplify = T)[,8], "_", simplify = T)[,1], 1, 1) %>% toupper()
    tmp = read_tsv(i) %>%
        mutate(type = paste(con, paste0("Mod", var), sep = "-")) %>%
        head(10)
    assign(paste(con, paste0("Mod", var),"GO", sep = "_"), tmp)
    print(dim(tmp))
    all_mod_go = rbind(all_mod_go, tmp)
}

all_mod_go = all_mod_go %>%
    mutate(EF = DOSE::parse_ratio(GeneRatio))
all_mod_go_filter = all_mod_go %>% filter(EF > 0.05)

all_mod_do_dup = table(all_mod_go_filter$ID) %>% as.data.frame() %>%
    `colnames<-`(c("ID", "freq")) %>%
    filter(freq > 1)


all_mod_go_filter %>%
    mutate(EF = DOSE::parse_ratio(GeneRatio),
           share = ifelse(ID %in% all_mod_do_dup$ID, "yes", "no")) %>% 
    arrange(desc(share)) %>%
    filter(pvalue < 0.005) %>%
    mutate(Description = factor(Description, levels = unique(Description))) %>%
    ggplot(aes(x = pvalue, y = Description, shape = share)) +
    geom_point(aes(size = EF, color = type, fill = type)) +
    scale_x_continuous(name = "p-value") +
    scale_y_discrete(name = "", labels = \(x) {str_wrap(x, width = 30)}) +
    scale_color_manual(values = col_nejm) +
    scale_fill_manual(values = alpha(col_nejm, 0.5)) +
    theme_minimal() +
    mytheme +
    theme(legend.position = "right",
          legend.background = element_blank(),
          axis.text.y = element_text(size = 8))

#ggsave(paste0(figdir, "FigS3/all_mod_go.pdf"), width = 180, height = 220, unit = "mm")
ggsave(paste0(figdir, "Figure6/all_mod_go_filter.pdf"), width = 130, height = 140, unit = "mm")

#all_mod_go %>% write_csv("../../../supplementary/20230210/SupplementaryTable2_plus.csv")

