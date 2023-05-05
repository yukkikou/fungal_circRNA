#!/usr/bin/Rscript

# signature
## Author: yuki
## Date: 2023.02.13

## environment
rm(list = ls())
setwd("D://_CircRNA/results/10_revise/20230201/")
rltdir="D://_CircRNA/results/"
set.seed(24)

## packages
source("D://_Scripts/R/_LIBRARY.R")
library(paletteer)
library(rlang)
library(ggpubr)

# colors and themes
source("D://_Scripts/R/_COLOR/ggsci_color_theme.R")

# functions
motifExact = function(flk_up, flk_down, protein_name){
    tmp_up = flk_up %>% filter(motif_struc %in% get(paste0(protein_name,"_motif"))) %>%
        select(circ_id) %>%
        mutate(flk_type = "up",
               protein = protein_name)
    tmp_down = flk_down %>% filter(motif_struc %in% get(paste0(protein_name,"_motif"))) %>%
        select(circ_id) %>%
        mutate(flk_type = "down",
               protein = protein_name)
    tmp = rbind(tmp_up, tmp_down)
    tmp$circ_id = gsub("\\(\\+\\)", "", tmp$circ_id)
    tmp
}

deExact = function(de, a){
    if (a == "all"){
        de %>% filter(pvalue < 0.05, abs(log2FoldChange) > 1) %>%
            select(feat_id) %>% unlist() %>% `names<-`(NULL)
    } else if (a == "up") {
        de %>% filter(pvalue < 0.05, log2FoldChange > 1) %>%
            select(feat_id) %>% unlist() %>% `names<-`(NULL)
    } else if (a == "down"){
        de %>% filter(pvalue < 0.05, log2FoldChange < (-1)) %>%
            select(feat_id) %>% unlist() %>% `names<-`(NULL)
    } else {
        stop("Input error.")
    }
}
# read in
circ_info = read_tsv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circRNA_info_all.tsv"))
circ_info_bed = read_tsv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circRNA_info_all.bed"), 
                         col_names = c("chr", "ss", "es", "circ_id_vertical", "len", "strand"))
circ_bsj = read_csv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circRNA_bsj.csv"))
circ_bsj_long = read_csv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circRNA_bsj_longer.csv"))
circ_ratio = read_csv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circRNA_ratio.csv"))
circ_ratio_long = read_csv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circRNA_ratio_longer.csv"))
circ_cpg = read_tsv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circ_CPG.tsv"))
circ_cpg_long = read_tsv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circ_CPG_longer.tsv"))
circ_de = read_csv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_circRNA_de.tsv")) %>%
    rename(circ_id_vertical = ...1)

circ_ratio_de = read_tsv(paste0(rltdir, "8_circ_WT_merge/7_expression/3_ratio/ratio/ratio_res.tsv")) %>%
    left_join(circ_info[,c(1,14)])

gene_exp = read_tsv(paste0(rltdir, "8_circ_WT_merge/merge/WT_merge_gene_circ_CPG.csv"))
gene_de = read_tsv(paste0(rltdir, "8_circ_WT_merge/7_expression/1_linear/gene/gene_res.tsv"))
pm1KnownGene = read_tsv("D://_CircRNA/scripts/Enrichment/gene.txt")

circ_splice = read_tsv(paste0(rltdir, "8_circ_WT_merge/4_spliceSite/all.motifRev.txt"))
circ_marker = read_tsv(paste0(rltdir, "8_circ_WT_merge/5_structure/circ_marker_table_type.tsv"))
circ_allinone = read_tsv(paste0(rltdir, "8_circ_WT_merge/5_structure/circ_allinone.tsv"))
circ_flk_up = read_tsv(paste0(rltdir, "8_circ_WT_merge/6_flanking/3_motif/1_streme/flk.up/sequences.tsv")) %>%
    filter(`motif_P-value` < 0.05, grepl("pm1_", seq_ID)) %>%
    `colnames<-`(c("motif_struc", "motif_id", "p", "circ_id", "seq_score", "seq_class", "holdout"))
circ_flk_down = read_tsv(paste0(rltdir, "8_circ_WT_merge/6_flanking/3_motif/1_streme/flk.down/sequences.tsv")) %>%
    filter(`motif_P-value` < 0.05, grepl("pm1_", seq_ID)) %>%
    `colnames<-`(c("motif_struc", "motif_id", "p", "circ_id", "seq_score", "seq_class", "holdout"))

motif_bind = read_tsv(paste0(rltdir, "8_circ_WT_merge/6_flanking/3_motif/2_tomtom/motif_enrich.tsv"))

# Section 1
# merge advances
bed_col = c("chr", "ss", "es", "circ_id", "score", "strand")
ciri2 = read_tsv(paste0(rltdir, "8_circ_WT_merge/1_number/CIRI2forCIRIquant.txt"), col_names = bed_col)
cirilong = read_tsv(paste0(rltdir, "8_circ_WT_merge/1_number/CIRIlongforCIRIquant.txt"), col_names = bed_col)
origin_intersect_circ = intersect(ciri2$circ_id, cirilong$circ_id) #784
merge_intersect_circ = circ_info$circ_id_vertical #3891
new_iden_circ = merge_intersect_circ[!(merge_intersect_circ %in% origin_intersect_circ)] #3132

new_iden_circ_cirilong = new_iden_circ[new_iden_circ %in% cirilong$circ_id] |> sort()
grep(new_iden_circ_cirilong[3], cirilong$circ_id, value = T)


sample_circ_bed = circ_info_bed[sample_circ_idx, ]
sample_circ_info = circ_info %>% filter(circ_id_vertical %in% sample_circ_bed$circ_id_vertical)
write_tsv(sample_circ_bed, "1_sampleCirc/sample_circ_15.bed", col_names = F)
write_tsv(sample_circ_info, "1_sampleCirc/sample_circ_15.info")
# reduce 15 to 10
#sam_10 = sample_circ_info[sample(1:15, size = 10, replace = F),]
#write_tsv(sam_10, "1_sampleCirc/sample_circ_15_to_10.info")
sam_10$circ_id_vertical

# output
circ_info %>% arrange(len) %>%
    head(5) %>%
    write_tsv("2_lenValidation/shortest_5.info")

circ_info %>% arrange(-len) %>%
    head(5) %>%
    write_tsv("2_lenValidation/longest_5.info")

circ_info_bed %>% filter(circ_id_vertical %in% select(head(arrange(circ_info, len), 5), circ_id_vertical)$circ_id_vertical) %>%
    write_tsv("2_lenValidation/shortest_5.bed", col_names = F)

circ_info_bed %>% filter(circ_id_vertical %in% select(head(arrange(circ_info, -len), 5), circ_id_vertical)$circ_id_vertical) %>%
    write_tsv("2_lenValidation/longest_5.bed", col_names = F)

# Section2
# linear isoform and circular isoform
gene_transnumber = read_tsv("D://_IGV/Talaromyces_marneffei/Annotation/gene_linear.number", col_names = c("linear_number", "gene_id"))
gene_circnumber = table(circ_info$gene_id) %>% as.data.frame() %>% `colnames<-`(c("gene_id", "circular_number"))
gene_isonumber = inner_join(gene_transnumber, gene_circnumber) %>% select(gene_id, linear_number, circular_number) %>%
    mutate(all_number = linear_number + circular_number,
           circular_per = round(circular_number/all_number*100,1))

# correlation
cor.test(gene_isonumber$circular_number, gene_isonumber$linear_number, method = "spearman")
# precentage quantiles
plot(x = seq(0, 1 ,0.05), y = quantile(gene_isonumber$circular_per, probs = seq(0,1,0.05)))
median(gene_isonumber$circular_per)    
mean(gene_isonumber$circular_per)    

# relationship
ggplot(gene_isonumber, aes(y = circular_per, x = linear_number)) +
    geom_point() +
    geom_smooth() +
    stat_regline_equation(label.y = 90, aes(label = ..eq.label..)) +
    stat_regline_equation(label.y = 80, aes(label = ..rr.label..))

ggplot(gene_isonumber, aes(y = circular_per, x = circular_number)) +
    geom_point() +
    geom_smooth() +
    stat_regline_equation(label.y = 90, aes(label = ..eq.label..)) +
    stat_regline_equation(label.y = 80, aes(label = ..rr.label..))
    
ggplot(gene_isonumber, aes(x = circular_number, y = linear_number)) +
    geom_point() +
    geom_smooth() +
    stat_regline_equation(label.y = 90, aes(label = ..eq.label..)) +
    stat_regline_equation(label.y = 80, aes(label = ..rr.label..))

# CPM
mean_circ_cpm = circ_cpg_long %>%
    group_by(circ_id,condition) %>%
    summarise(mean = mean(cpg)/1000) %>%
    pivot_wider(names_from = condition, 
                names_prefix = "cpm_",
                values_from = mean) %>%
    mutate(low_mycelium = cpm_mycelium < 0.01,
           low_yeast = cpm_yeast < 0.01,
           low_both = low_mycelium & low_yeast,
           sum_mean = cpm_mycelium + cpm_yeast) %>%
    arrange(-sum_mean)
# ratio
mean_circ_ratio = circ_ratio_long %>%
    group_by(circ_id,condition) %>%
    summarise(mean = mean(ratio)) %>%
    pivot_wider(names_from = condition, 
                names_prefix = "ratio_",
                values_from = mean) %>%
    mutate(low_mycelium = ratio_Mycelium < 0.1,
           low_yeast = ratio_Yeast < 0.1,
           low_both = low_mycelium & low_yeast,
           sum_mean = ratio_Mycelium + ratio_Yeast) %>%
    arrange(-sum_mean)

# low expression circRNA
length(intersect(mean_circ_ratio[mean_circ_ratio$low_both, 'circ_id']$circ_id,
          mean_circ_cpm[mean_circ_cpm$low_both, 'circ_id']$circ_id))

# high expression circRNAs
intersect(mean_circ_cpm[mean_circ_cpm$sum_mean>0.5,]$circ_id,
          mean_circ_ratio[mean_circ_ratio$sum_mean<2 & mean_circ_ratio$sum_mean>1,]$circ_id)

circ_marker %>% filter(circ_id_vertical %in% c("pm1_04:1410708|1410903","pm1_02:869646|869895")) %>%
    select(circ_id_vertical, circ_label, gene_id, circ_marker)

mean_circ_cpm %>% filter(circ_id %in% c("pm1_04:1410708|1410903","pm1_02:869646|869895"))
mean_circ_ratio %>% filter(circ_id %in% c("pm1_04:1410708|1410903","pm1_02:869646|869895"))


circ_info %>% filter(circ_id_vertical == "pm1_01:2332393|2332684") %>% # TM010856.00: EF-1 alpha subunit
    select(circ_id_vertical, circ_label, gene_id)

mean_circ_cpm %>% filter(circ_id == "pm1_01:2332393|2332684")
# Section3


# Section4
rr_circ= "pm1_02:2839992|2840162"
circ_info %>% filter(circ_id_vertical == rr_circ) %>% as.data.frame()

# Section5
hnRNP_motif = unique(unlist(motif_bind[grepl("HNRNP", motif_bind$gene_id),"motif"]))
SF_motif = unique(unlist(motif_bind[grepl("SF", motif_bind$gene_id),"motif"]))

hnRNP_circ = motifExact(circ_flk_up, circ_flk_down, "hnRNP")
SF_circ = motifExact(circ_flk_up, circ_flk_down, "SF")
circ_id_overlap = intersect(hnRNP_circ$circ_id, SF_circ$circ_id)

circ_rnp_overlap = rbind(hnRNP_circ, SF_circ)
circ_rnp_overlap$protein = ifelse(circ_rnp_overlap$circ_id %in% circ_id_overlap, 
                                  "both", circ_rnp_overlap$protein)

table(circ_rnp_overlap$protein, circ_rnp_overlap$flk_type)

# Section6
gene_de_label = gene_de %>% filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    mutate(de_gene_type = ifelse(log2FoldChange > 1, "ug", "dg")) %>%
    right_join(gene_de) %>%
    mutate(de_gene_type = ifelse(is.na(de_gene_type), "ng", de_gene_type))

circ_de_label = circ_de %>% left_join(circ_info) %>% 
    filter(DE != 0) %>%
    mutate(de_circ_type = ifelse(logFC > 1, "uc", "dc")) %>%
    select(circ_id_vertical:len, gene_id, circtrans_id, de_circ_type)
    
circ_ratio_de_label = circ_ratio_de %>%
    filter(pvalue < 0.05, abs(log2FoldChange) > 1) %>%
    mutate(de_ratio_type = ifelse(log2FoldChange > 1, "ur", "dr"))

circ_de_dt = inner_join(circ_de_label, circ_ratio_de_label) %>%
    inner_join(gene_de_label, by = c("gene_id" = "gene_id"), suffix = c(".ratio", ".gene")) %>%
    mutate(de_cls_type = paste(de_ratio_type, de_gene_type, de_circ_type, sep = "-"))


circ_dd_sig_label = data.frame(
    "classification" = c(rep("intensified", 4), rep("buffered", 2), rep("reverse", 3)),
    "de_cls_type" = c("dr-ng-dc", "ur-ug-uc", "ur-dg-dc", "ur-ng-uc",
                      "dr-ug-uc", "dr-dg-dc", 
                      "dr-dg-uc", "dr-ug-dc", "ur-dg-uc")
)
circ_clasf_map = data.frame(
    "de_cls_type" = paste(paste(rep(paste0(c("u", "d", "n"), "r"), 3),
                                rep(paste0(c("u", "d", "n"), "g"), each = 3), sep = "-"),
                          rep(paste0(c("u", "d", "n"), "c"), each = 9), sep = "-") %>% sort()
) %>% left_join(circ_dd_sig_label) %>%
    mutate(classification = ifelse(is.na(classification), "background", classification)) %>%
    right_join(circ_de_dt) %>%
    #left_join(pm1KnownGene) %>%
    as_tibble()

table(circ_clasf_map$classification)
circ_clasf_map %>% filter(!is.na(gene_name)) %>% as.data.frame()
circ_clasf_map %>% write_tsv("3_de_reAnalysis/circ_clasf_map.tsv")
circ_clasf_map %>% write_tsv(paste0(rltdir, "8_circ_WT_merge/8_DE/circ_clasf_map.tsv"))

# rank
View(circ_clasf_map)
circ_clasf_map %>% filter(grepl("1328636", circ_id_vertical)) %>% as.data.frame()
final_de_circ_table = circ_clasf_map %>% filter(!is.na(len)) %>%
    left_join(mean_circ_cpm, by = c("circ_id_vertical" = "circ_id")) %>%
    left_join(mean_circ_ratio, by = c("circ_id_vertical" = "circ_id"), suffix = c(".cpg", ".ratio")) %>%
    filter(classification != "buffered") %>%
    filter(!low_both.ratio, !low_both.cpg) %>%
    filter(!low_yeast.ratio, !low_mycelium.ratio) %>%
    filter(!(sum_mean.ratio %in% 1:2)) %>%
    mutate(score = FDR/(sum_mean.cpg*sum_mean.ratio*2^logFC*2^log2FoldChange.ratio)) %>%
    arrange(score)

final_de_circ_table %>%
    write_tsv(paste0(rltdir, "8_circ_WT_merge/8_DE/final_de_circ_table.tsv"))

final_de_circ_table %>%
    select(circ_id_vertical) %>%
    left_join(circ_info) %>%
    write_tsv(paste0(rltdir, "8_circ_WT_merge/8_DE/final_de_circ_table.info"))

final_de_circ_table %>%
    select(circ_id_vertical) %>%
    left_join(circ_info_bed) %>%
    write_tsv(paste0(rltdir, "8_circ_WT_merge/8_DE/final_de_circ_table.bed"), col_names = F)


# Section7
# analysis in results/network

dedir="D:/_CircRNA/results/9_WGCNA/1_RRHO/1_DEseq2/"
mkw_de = read_tsv(paste0(dedir, "mycelium_kd_wt_res.tsv"))
ykw_de = read_tsv(paste0(dedir, "yeast_kd_wt_res.tsv"))
mow_de = read_tsv(paste0(dedir, "mycelium_oe_wt_res.tsv"))
yow_de = read_tsv(paste0(dedir, "yeast_oe_wt_res.tsv"))

# mark de_ud_gene and de_ud_circRNA number
# include overlap (ud)
mudir = "D://_CircRNA/results/7_mutants/"
geneMtx = read_tsv(paste0(mudir, "2_linear/linear_gene_count.clean.tsv"))
circMtx = read_csv(paste0(mudir, "1_circ/merge_circRNA_bsj.csv")) %>%
    dplyr::rename(feat_id = circ_id)

for (i in c("m", "y")){
    for (j in c("k", "o")){
        tmp = get(paste0(i, j, "w_de"))
        tmp = tmp %>% mutate(de_type = ifelse(feat_id %in% deExact(tmp, a = "all"), 
                                           ifelse(feat_id %in% deExact(tmp, a = "up"), "up", "down"), NA))
        assign(paste0(i, j, "w_de"), tmp)
        print(table(get(paste0(i, j, "w_de"))$de_type))
    }
}

mu_de_table = data.frame(
    "feat_id" = c(geneMtx$Geneid, circMtx$feat_id)
) %>% full_join(
    full_join(full_join(mkw_de[,c(1,3,6:8)], ykw_de[,c(1,3,6:8)], 
                                    by = c("feat_id" = "feat_id"),
                                    suffix = c(".mkw", ".ykw")),
                          full_join(mow_de[,c(1,3,6:8)], yow_de[,c(1,3,6:8)], 
                                    by = c("feat_id" = "feat_id"),
                                    suffix = c(".mow", ".yow")))
) %>% as_tibble()

tmp_df = data.frame(matrix(nrow = 0, ncol = 6))

for (i in c("m", "y")){
    for (j in c("k", "o")){
        for (k in c("up", "down")){
            var = paste0(i, j, "w")
            tmp_de = mu_de_table[mu_de_table[paste0("de_type.", var)] == k, ] %>% filter(!is.na(feat_id))
            
            tmp = mu_de_table[mu_de_table[paste0("de_type.", var)] == k, ] %>%
                filter(!is.na(feat_id)) %>%
                select(feat_id, contains("de_type"), -contains(var))
            tmp_gene = tmp %>% filter(grepl("TM", feat_id))
            tmp_circ = tmp %>% filter(grepl("pm", feat_id))
            
            n_gene = nrow(tmp_gene[apply(tmp_gene[,-1], 1, function(y) any(!is.na(y))),])
            n_circ = nrow(tmp_circ[apply(tmp_circ[,-1], 1, function(y) any(!is.na(y))),])
            
            tmp_k = data.frame("mutant" = var, type = k)
            tmp_k[,"de_gene"] = length(grep("TM", tmp_de$feat_id))
            tmp_k[,"de_circ"] = length(grep("pm", tmp_de$feat_id))
            tmp_k[,"over_gene"] = n_gene
            tmp_k[,"over_circ"] = n_circ
            tmp_df = rbind(tmp_df, tmp_k)
        }
    }
}


mu_df = tmp_df %>% mutate(
    nonover_gene = de_gene - over_gene,
    nonover_circ = de_circ - over_circ
) %>% pivot_longer(cols = de_gene:nonover_circ,
                        names_to = "number_type",
                        values_to = "number") %>%
    filter(!grepl("de", number_type)) %>%
    mutate(feat = str_split(number_type, "_", simplify = T)[,2],
           trend = str_split(number_type, "_", simplify = T)[,1],
           type = paste(type, feat, sep = "_")) %>%
    select(mutant, type, trend, number)

# WGCNA network
head(pm1KnownGene)
clist = list.files(path = wcgnsdir, pattern = "CytoscapeInput-edges-*", recursive = T, full.names = T)
all_mod_cyto = data.frame(matrix(ncol = 7, nrow = 0)) 
for (i in clist){
    var = substr(str_split(str_split(i, "/", simplify = T)[,8], "-", simplify = T)[,3], 1, 1) %>% toupper()
    con = substr(str_split(i, "/", simplify = T)[,7], 1, 1) %>% toupper()
    tmp = read_tsv(i) %>%
        mutate(type = paste(con, paste0("Mod", var), sep = "-")) %>%
        arrange(desc(weight)) %>%
        head(30)
    assign(paste(con, paste0("Mod", var),"Cyto", sep = "_"), tmp)
    print(dim(tmp))
    all_mod_cyto = rbind(all_mod_cyto, tmp)    
}

all_mod_cyto = all_mod_cyto %>%
    left_join(pm1KnownGene, by = c("fromNode" = "gene_id")) %>%
    left_join(pm1KnownGene, by = c("toNode" = "gene_id")) %>%
    mutate(fromAltName = gene_name.x,
           toAltName = gene_name.y) %>%
    select(-contains("gene_name"))

for (i in unique(all_mod_cyto$type)){
    all_mod_cyto %>%
        filter(type == i) %>%
        write_tsv(paste0(wcgnsdir, "/top30forplot/", i, ".txt"))
}

# circDS-1 in mycelium
#all_geneMod %>% filter(grepl("1328808",feat_id))
M_circGene = read_tsv(paste0(wcgnsdir, "/circDS-1/Cyto-edges-Mgrey60_circDS.txt")) %>% 
    arrange( desc(weight)) %>% 
    filter((grepl("1328808",fromNode)|grepl("1328808",toNode))) 
#%>%    filter(weight >= 0.6)

M_circEdges = read_tsv(paste0(wcgnsdir, "/circDS-1/Cyto-edges-Mgrey60_circDS.txt")) %>% 
    arrange(desc(weight))

for (i in seq_along(M_circGene$toNode)){
    M_circGene = rbind(M_circGene, cytoExt10(M_circEdges, M_circGene$toNode[i]))
}

M_circGene = M_circGene[!duplicated(M_circGene),]
write_tsv(M_circGene, paste0(wcgnsdir, "/circDS-1/m-circ-10.tsv"))


Y_circGene = read_tsv(paste0(wcgnsdir, "/circDS-1/Cyto-edges-Ylightgreen_circDS.txt")) %>% 
    arrange( desc(weight)) %>% 
    filter((grepl("1328808",fromNode)|grepl("1328808",toNode))) 
#%>% filter(weight >= 0.3)

Y_circEdges = read_tsv(paste0(wcgnsdir, "/circDS-1/Cyto-edges-Ylightgreen_circDS.txt")) %>% 
    arrange(desc(weight))


for (i in seq_along(Y_circGene$toNode)){
    Y_circGene = rbind(Y_circGene, cytoExt10(Y_circEdges, Y_circGene$toNode[i]))
}
Y_circGene= Y_circGene[!duplicated(Y_circGene),]
write_tsv(Y_circGene, paste0(wcgnsdir, "/circDS-1/Y-circ-10.tsv"))

# gene annotation
pm1Omics |> head()

M_circGene_anno = M_circGene %>% left_join(pm1Omics, by = c("fromNode" = "SequenceName")) %>%
    left_join(pm1Omics, by = c("toNode" = "SequenceName")) %>%
    mutate(fromAltName = SequenceDescription.x,
           toAltName = SequenceDescription.y,
           fromAltName = ifelse(is.na(fromAltName), fromNode, fromAltName),
           toAltName = ifelse(is.na(toAltName), toNode, toAltName)) %>%
    select(-contains("Sequence"))

write_tsv(M_circGene_anno, paste0(wcgnsdir, "/circDS-1/m-circ-all-anno.tsv"))

Y_circGene_anno = Y_circGene %>% left_join(pm1Omics, by = c("fromNode" = "SequenceName")) %>%
    left_join(pm1Omics, by = c("toNode" = "SequenceName")) %>%
    mutate(fromAltName = SequenceDescription.x,
           toAltName = SequenceDescription.y,
           fromAltName = ifelse(is.na(fromAltName), fromNode, fromAltName),
           toAltName = ifelse(is.na(toAltName), toNode, toAltName)) %>%
    select(-contains("Sequence"))

write_tsv(Y_circGene_anno, paste0(wcgnsdir, "/circDS-1/y-circ-all-anno.tsv"))