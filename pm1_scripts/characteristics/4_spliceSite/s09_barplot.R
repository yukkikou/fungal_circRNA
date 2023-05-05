#!/usr/bin/Rscript

setwd("/media/data5/hxy/circResult/4_spliceSite")

library(tidyverse)
library(RColorbrewer)
library(patchwork)

cols = brewer.pal(8, "Set2")
mytheme = theme(
	plot.title = element_text(hjust = 0.5),
	text = element_text(size = 10),
	axis.text.x = element_text(size = 8, color = "black", hjust = 0.5, vjust = 0.5, angle = 45),
	axis.text.y = element_text(size = 8, color = "black"),
	legend.title = element_text(size = 8, color = "grey33"),
	legend.text = element_text(size = 6, color = "grey33"),
	legend.key.size = unit(11, "pt"),
	legend.position = "bottom",
	panel.grid.major=element_blank(),panel.grid.minor=element_blank()
)

circ_ss_table = read_tsv("all.motifRev.table", col_names = c("motif", "count", "freq")) %>%
	mutate(type = "circRNA") %>%
  filter(count > 1)
exon_ss_table = read_tsv("pm1_exon/pm1_intron.motifRev.table", col_names = c("motif", "count", "freq")) %>%
	mutate(type = "intron") %>%
  filter(count > 1)


plot_circ_motif_table_compare_non = ggplot() +
    geom_hline(yintercept = 0, color = "grey60", size = 0.3) +
    geom_col(data = filter(circ_ss_table, motif != "AGGT"),
        mapping = aes(x = motif, y = round(freq*100,2), fill = type),
        width = 0.4) +
    geom_text(data = filter(circ_ss_table, motif != "AGGT"),
        mapping = aes(x = motif, y = round(freq*100,2), label = round(freq*100,2)),
        size = 2, color = "grey60",
        vjust = -1) +
    geom_col(data = filter(exon_ss_table, motif != "AGGT"),
        mapping = aes(x = motif, y = -round(freq*100,2), fill = type),
        width = 0.4) +
    geom_text(data = filter(exon_ss_table, motif != "AGGT"),
        mapping = aes(x = motif, y = -round(freq*100,2), label = round(freq*100,2)),
        size = 2, color = "grey60",
        vjust = 1.5) +
    scale_fill_manual(values = cols) +
    labs(x = "", y = "Precentage of motifs") +
    scale_y_continuous(limits = c(-2.5,3.5), breaks = c(-2,-1,0,1,2), labels = c(2,1,0,1,2)) +
    theme_bw()+
    mytheme

plot_circ_motif_table_compare_can = ggplot() +
  geom_hline(yintercept = 0, color = "grey60", size = 0.3) +
  geom_col(data = filter(circ_ss_table, motif != "AGGT"),
           mapping = aes(x = motif, y = round(freq*100,2), fill = type),
           width = 0.4) +
  geom_text(data = filter(circ_ss_table, motif != "AGGT"),
            mapping = aes(x = motif, y = round(freq*100,2), label = round(freq*100,2)),
            size = 2, color = "grey60",
            vjust = -1) +
  geom_col(data = filter(exon_ss_table, motif != "AGGT"),
           mapping = aes(x = motif, y = -round(freq*100,2), fill = type),
           width = 0.4) +
  geom_text(data = filter(exon_ss_table, motif != "AGGT"),
            mapping = aes(x = motif, y = -round(freq*100,2), label = round(freq*100,2)),
            size = 2, color = "grey60",
            vjust = 1.5) +
  scale_fill_manual(values = cols) +
  labs(x = "", y = "Precentage of motifs") +
  scale_y_continuous(limits = c(-2.5,3.5), breaks = c(-2,-1,0,1,2), labels = c(2,1,0,1,2)) +
  theme_bw()+
  mytheme


ggsave("circ_motif_table_compare.pdf", width = 150, height = 120, units = "mm")


