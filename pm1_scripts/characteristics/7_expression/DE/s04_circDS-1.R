setwd("D://_CircRNA/results/5_circRNA_merging/7_expression/")

circ_table |> head()

de_circ_res |> head()
de_ratio_res |> head()

circ_bsj |> head()
circ_cpm |> head()

circ_full_length = read_tsv("../2_annotation/circ_IN.length", col_names = c("circ_id", "type", "length"))

# full length selection
circ_selection = circ_table %>% filter(circ_id %in% circ_full_length$circ_id) %>%
  left_join(circ_full_length, by = c("circ_id" = "circ_id",
                                     "circ_type" = "type"))

# CPM selection
circ_cpm_sum = data.frame(
  mycelium_mean_cpm = apply(circ_cpm, 1, \(x){mean(x[1:4])}),
  yeast_mean_cpm = apply(circ_cpm, 1, \(x){mean(x[5:8])}),
  all_mean_cpm = apply(circ_cpm, 1, \(x){mean(x)})
)

circ_cpm_sum$max = apply(circ_cpm_sum, 1, \(x){max(x)})
circ_cpm_sum = circ_cpm_sum %>%
  rownames_to_column("circ_id") %>%
  arrange(-max)

circ_selection = circ_selection %>% left_join(circ_cpm_sum) %>%
  arrange(-max)

# add circRNA de information
circ_selection = circ_selection %>% left_join(de_circ_res) %>%
  filter(DE != 0) 

# add ratio de selection
circ_selection = circ_selection %>% left_join(de_ratio_res[,c(1,3,7)])

# add ratio
circ_ratio  = read_csv("IN_merge_circRNA_ratio.csv") %>%
  column_to_rownames("circ_id")
circ_ratio
circ_ratio_sum = data.frame(
  mycelium_mean_ratio = apply(circ_ratio, 1, \(x){mean(x[1:4])}),
  yeast_mean_ratio = apply(circ_ratio, 1, \(x){mean(x[5:8])}),
  all_mean_ratio = apply(circ_ratio, 1, \(x){mean(x)})
)

circ_ratio_sum$max = apply(circ_ratio_sum, 1, \(x){max(x)})
circ_ratio_sum = circ_ratio_sum %>%
  rownames_to_column("circ_id") %>%
  arrange(-max)

circ_selection = circ_selection %>% left_join(circ_ratio_sum, 
                                              by = c("circ_id" = "circ_id"),
                                              suffix = c(".cpm", ".ratio"))

# final rank
circ_selection_final = circ_selection %>%
  filter(FDR<0.05, padj<0.05,abs(logFC)>1, abs(log2FoldChange)>1,max.ratio>0.5, max.ratio<1) %>%
  arrange(-abs(logFC),FDR,  -max.ratio,-abs(log2FoldChange),  padj) %>%
  select(circ_id, de_cls_type, length, mycelium_mean_ratio, yeast_mean_ratio, max.ratio, padj, FDR, logFC, log2FoldChange, max.cpm)

# add 
colnames(circ_selection_final)

circ_selection_final = circ_selection_final %>% left_join(circ_to_gene)
table(circ_selection_final$gene_id %in% pm1ReportedGene)

circ_selection_final %>% left_join(pm1GeneGo[!duplicated(pm1GeneGo[,c(1,3,4,5)]), c(1,3,4,5)], by = c("gene_id" = "Gene")) %>% 
  write_tsv("circ_selection_final_list.tsv")


