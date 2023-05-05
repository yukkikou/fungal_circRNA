#!/usr/bin/Rscript
setwd("//211.82.71.164/media/data5/hxy/circResult/2_annotation")

library(tidyverse)
library(RColorBrewer)
library(ggbreak)

cols = brewer.pal(8, "Set2")
mytheme = theme(plot.title = element_text(hjust = 0.5),
                text=element_text(size=10),
                axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
                axis.text.y = element_text(size = 8,color = "black"),
                legend.text = element_text(size=6,color = "black"),
                legend.title = element_text(size=8),
                legend.key.size = unit(11, "pt"))

# read in
circ_gene_map = read_tsv("circ_gene.map", col_names = c("circ_id", "gene_id"))

myc_circ_id = readLines("../1_number/IN_M.id")
yst_circ_id = readLines("../1_number/IN_Y.id")
di_circ_id = intersect(myc_circ_id, yst_circ_id)

lth = read_tsv("circ_IN.length", col_names = c("circ_id", "type", "length")) %>%
  mutate(source = NA)

for(i in 1:dim(lth)[1]){
  lth[i,"source"] = ifelse(lth[i,"circ_id"] %in% di_circ_id, "Both",
                           ifelse(lth[i,"circ_id"] %in% myc_circ_id, "Mycelium","Yeast"))
}

# visualization
plot_circ_gene_iso = circ_gene_map$gene_id %>% table() %>%
  table() %>% as.data.frame() %>%
  `colnames<-`(c("iso","freq")) %>%
  ggplot(aes(x = iso, y = freq, label = freq)) +
    geom_col(width = 0.4, fill = cols[1]) +
    geom_text(size = 3, color = "grey33", hjust = -0.25) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1500)) +
    labs(x = "CircRNA isoform number (per gene)", y = "Frequency") +
    coord_flip() +
    theme_bw() +
    mytheme
ggsave("circ_gene_iso.pdf", width = 70, height = 120, units = "mm")

plot_circ_full_length_source = lth %>% 
  ggplot(aes(x = length, color = source)) +
		geom_density() +
    scale_color_manual(values = cols) +
    labs(x = "Length of circRNA", y = "Density") +
		theme_bw() +
    mytheme +
    theme(legend.position = "top")
ggsave("circ_full_length.pdf", width = 70, height = 70, units = "mm")

plot_circ_full_length_type = lth %>% 
  ggplot(aes(x = length, color = type)) +
  geom_density() +
  facet_wrap(~source) +
  scale_color_manual(values = cols) +
  labs(x = "Length of circRNA", y = "Density") +
  theme_bw() +
  mytheme +
  theme(legend.position = "top")
ggsave("circ_full_length_type.pdf", width = 140, height = 70, units = "mm")
