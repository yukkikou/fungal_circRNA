source("D://_CircRNA/scripts/Enrichment/s00_goFunction.R")

# input
gene_de = read_tsv("../../5_circRNA_merging/7_expression/_linear_by_featureCounts/both/gene/gene_res_sig.tsv")
gene_de_up = read_tsv("../../5_circRNA_merging/7_expression/_linear_by_featureCounts/both/gene/gene_res_up.tsv")
gene_de_down = read_tsv("../../5_circRNA_merging/7_expression/_linear_by_featureCounts/both/gene/gene_res_down.tsv")

circ_de = read_csv("../../5_circRNA_merging/7_expression/IN_merge_circRNA_de.tsv")
colnames(circ_de)[1] = "circ_id"
circ_info = read_csv("../../5_circRNA_merging/7_expression/IN_merge_circRNA_info.csv")

circ_map = read_tsv("../../5_circRNA_merging/7_expression/ratio/circ_to_gene.tsv")
ratio_de = read_tsv("../../5_circRNA_merging/7_expression/ratio/ratio/ratio_res_sig.tsv") %>% 
  left_join(circ_map, by = c("circtrans_id" = "gene_name"))
ratio_de_up = read_tsv("../../5_circRNA_merging/7_expression/ratio/ratio/ratio_res_up.tsv") %>% 
  left_join(circ_map, by = c("circtrans_id" = "gene_name"))
ratio_de_down = read_tsv("../../5_circRNA_merging/7_expression/ratio/ratio/ratio_res_down.tsv") %>% 
  left_join(circ_map, by = c("circtrans_id" = "gene_name"))


# de circRNA with gene id
circ_sig_all = circ_de %>% filter(DE != 0) %>%
  left_join(circ_info, by = c("circ_id" = "circ_id"))

circ_sig_up = circ_sig_all %>% filter(DE == 1)
circ_sig_down = circ_sig_all %>% filter(DE == -1)

# enrichment
circ_res = go_result(circ_sig_all$gene_id) %>% filter(p.adjust < 0.05)
circ_res_up = go_result(circ_sig_up$gene_id) %>% filter(p.adjust < 0.05)
circ_res_down = go_result(circ_sig_down$gene_id) %>% filter(p.adjust < 0.05)

gene_de_go = go_result(gene_de$gene_id) %>% filter(p.adjust < 0.05)
gene_de_up_go = go_result(gene_de_up$gene_id) %>% filter(p.adjust < 0.05)
gene_de_down_go = go_result(gene_de_down$gene_id) %>% filter(p.adjust < 0.05)


ratio_de_go = go_result(ratio_de$gene_id) %>% filter(pvalue < 0.01)
ratio_de_up_go = go_result(ratio_de_up$gene_id) %>% filter(pvalue < 0.01)
ratio_de_down_go = go_result(ratio_de_down$gene_id) %>% filter(pvalue < 0.01)

# output
write_tsv(circ_res, "../../results/5_circRNA_merging/7_expression/circ_de_go.tsv")
write_tsv(gene_de_go, "../../results/5_circRNA_merging/7_expression/gene_de_go.tsv")
write_tsv(gene_de_up_go, "../../results/5_circRNA_merging/7_expression/gene_de_up_go.tsv")
write_tsv(gene_de_down_go, "../../results/5_circRNA_merging/7_expression/gene_de_down_go.tsv")

# visualization
# go plot
#circ_res %>%
gene_de_up_go %>%
#gene_de_down_go %>%
  ggplot(aes(x = Description, y = DOSE::parse_ratio(GeneRatio)))+
  geom_point(aes(size = Count, color=(-log10(p.adjust))))+
  coord_flip()+
  scale_x_discrete(labels=function(x) stringr::str_wrap(x, width=50))+
  scale_color_continuous(type = "viridis") +
  theme_bw()+
  labs(x=NULL,y="Gene ratio")+
  guides(color=guide_legend(title="-log10(FDR)"))+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=10),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        legend.text = element_text(size=6,color = "black"),
        legend.title = element_text(size=8),
        legend.key.size = unit(11, "pt"))

# volcano
circ_de %>% 
  ggplot(aes(x = logFC, y = -log10(FDR)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  geom_point(aes(size=-log10(FDR), color= -log10(FDR)))+
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

ratio_de_down_go %>%
  #gene_de_down_go %>%
  ggplot(aes(x = Description, y = DOSE::parse_ratio(GeneRatio)))+
  geom_point(aes(size = Count, color=(-log10(pvalue))))+
  coord_flip()+
  scale_x_discrete(labels=function(x) stringr::str_wrap(x, width=50))+
  scale_color_continuous(type = "viridis") +
  theme_bw()+
  labs(x=NULL,y="Gene ratio")+
  guides(color=guide_legend(title="-log10(pvalue)"))+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=10),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        legend.text = element_text(size=6,color = "black"),
        legend.title = element_text(size=8),
        legend.key.size = unit(11, "pt"))

