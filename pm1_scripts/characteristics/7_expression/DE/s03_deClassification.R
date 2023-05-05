# classfication de circ/linear type due to ratio change
setwd("D://_CircRNA/results/5_circRNA_merging/7_expression/")
#source("../../../scripts/Enrichment/s00_goFunction.R")
source("D://_Scripts/R/_COLOR/colors.R")

library(tidyverse)
library(RColorBrewer)

cols = brewer.pal(12, "Set3")[-2]
mytheme = theme(plot.title = element_text(hjust = 0.5),
                text=element_text(size=10),
                axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
                axis.text.y = element_text(size = 8,color = "black"),
                legend.text = element_text(size=6,color = "black"),
                legend.title = element_text(size=8),
                legend.key.size = unit(11, "pt"))

# prior variables (do not change)
ls()
de_ratio_res |> head()
de_circ_res |> head()
de_gene_res |> head()
circ_to_gene |> head()

# new table ready to generated
circ_table = circ_info[,1:3] %>% as_tibble() %>%
  inner_join(circ_to_gene) %>% 
  dplyr::rename(circtrans_id = gene_name) %>%
  mutate(de_ratio_type = NA,
         de_circ_type = NA,
         de_gene_type = NA,
         de_cls_type = NA)
nrCirctable = nrow(circ_table)

# non_sig: adjp<0.05 & abs(log2FC)>1
non_ratio_circ = de_ratio_res %>% filter(padj >= 0.05| abs(log2FoldChange) <= 1)
non_gene_circ = de_gene_res %>% filter(padj >= 0.05 | abs(log2FoldChange) <= 1)
circ_res_down = de_circ_res %>% filter(DE == (-1))
circ_res_up = de_circ_res %>% filter(DE == 1)
circ_res_non = de_circ_res %>% filter(DE == 0)

for (i in 1:nrCirctable){
  circtrans_tmp = circ_table[i, "circtrans_id"]
  circ_table[i, "de_ratio_type"] = ifelse(circtrans_tmp %in% c(non_ratio_circ$circtrans_id,
                                                               ratio_res_na$circtrans_id), "nr", 
                                          ifelse(circtrans_tmp %in% ratio_res_up$circtrans_id,
                                                 "ur", "dr"))
  
  circid_tmp = circ_table[i, "circ_id"]
  circ_table[i, "de_circ_type"] = ifelse(circid_tmp %in% c(circ_res_non$circ_id), "nc", 
                                          ifelse(circid_tmp %in% circ_res_up$circ_id,
                                                 "uc", "dc"))
  geneid_tmp = circ_table[i, "gene_id"]
  circ_table[i, "de_gene_type"] = ifelse(geneid_tmp %in% c(non_gene_circ$gene_id,
                                                           gene_res_na$gene_id), "ng", 
                                         ifelse(geneid_tmp %in% gene_res_up$gene_id,
                                                "ug", ifelse(geneid_tmp %in% gene_res_down$gene_id,"dg",NA)))
  
  circ_table[i, "de_cls_type"] = paste(circ_table[i, "de_ratio_type"],
                                       circ_table[i, "de_gene_type"],
                                       circ_table[i, "de_circ_type"], sep = "-")
  
}

# statistics
table(circ_table$de_ratio_type)
table(circ_table$de_circ_type)
table(circ_table$de_gene_type)
table(circ_table$de_cls_type) |> as.data.frame() 

# subgroups
stable_circ_pairs = circ_table %>% filter(de_cls_type %in% c("nr-ug-uc","nr-dg-dc"))
go_result(stable_circ_pairs$gene_id) %>% filter(p.adjust<0.05) %>% select(Description)

intensifirf_circ_paird = circ_table %>% filter(de_cls_type %in% c("ur-ng-uc","ur-ug-uc", "dr-ng-dc","dr-dg-dc"))
go_result(intensifirf_circ_paird$gene_id) %>% filter(p.adjust<0.05) %>% select(Description)

reverse_circ_pairs = circ_table %>% filter(de_cls_type %in% c("ur-dg-nc","ur-dg-uc","ur-dg-dc",
                                                              "dr-ug-nc","dr-ug-uc","dr-ug-dc"))
go_result(reverse_circ_pairs$gene_id) %>% filter(p.adjust<0.05) %>% select(Description)

# reported gene overlap
pm1ReportedGene = read_tsv("../../../scripts/Enrichment/gene.txt")
circ_table %>% inner_join(pm1ReportedGene, by = c("gene_id" = "gene_id")) %>%
  arrange(de_cls_type) %>%
  write_tsv("circ_table_pm1ReportedGene.tsv")

# isoform bar
plot_ratio_iso = circ_table$gene_id %>% table() %>% table() %>%
  as.data.frame() %>% `colnames<-`(c("iso", "freq")) %>%
  ggplot(aes(x = iso, y = freq, label = freq)) +
    geom_col(position = "dodge", fill = cols[3], width = 0.6) +
    geom_text(vjust = -0.5, size = 2.5, color = "grey33") +
    labs(x = "CircRNA isoform number", y = "Frequency") +
    theme_bw() +
  mytheme

# volcano plot
de_gene_res_anno = de_gene_res %>% inner_join(circ_to_gene, by = c("gene_id" = "gene_id")) %>%
  select(gene_id, log2FoldChange, padj, circ_id, gene_name) %>%
  dplyr::rename(geneFC = log2FoldChange,
                geneFDR = padj)
de_circ_res_anno = de_circ_res %>% inner_join(circ_to_gene, by = c("circ_id" = "circ_id")) %>%
  select(circ_id, logFC, FDR, gene_id, gene_name) %>%
  dplyr::rename(circFC = logFC,
                circFDR = FDR)

de_gene_circ_table = left_join(de_gene_res_anno, de_circ_res_anno,
                               by = c("gene_name" = "gene_name",
                                      "gene_id" = "gene_id",
                                      "circ_id" = "circ_id")) %>%
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
    geom_point(size = 1.8, alpha = 0.65) +
    xlim(c(-11,11)) +
    ylim(c(-11,11)) +
    scale_color_manual(values = c(col8[c(5,6,7)], "grey60")) +
    labs(x = "Log2(FoldChange) of gene", y = "Log2(FoldChange) of circRNA") +
    theme_bw() +
    mytheme +
    theme(
      legend.position = "none",
      panel.grid=element_blank()
    )

save(list = ls(), "DE.RData")
