### PLOT SPLICE SITE OF CIRCRNA

setwd("//211.82.71.8/media/data7/hxy/PM1_illumima_202101/ANNO_XS01KF2020120029_PM-XS01KF2020120029-02_AHF3JTCCX2_2021-01-08/Rawdata/6_characteristics/2_splitSite/stat")

source("D://_Scripts/R/_COLOR/theme.R")

library(ggplot2)
library(RColorBrewer)

load("sum_stat.RData")
sum_stat = sum_stat[-1,]

sum_stat |> head()
lf = dplyr::filter(sum_stat, Sample == "All") |>
  dplyr::select(Count, Motif) |> 
  dplyr::arrange(Count) |> 
  dplyr::select(Motif) 

lab = data.frame(
  ID = c("All", paste0("IN_M",1:4), paste0("IN_Y",c(3:4,1:2))),
  Sample = c(sort(unique(sum_stat$Sample)))
)

sum_stat = dplyr::left_join(sum_stat, lab)

sum_stat$Motif = factor(sum_stat$Motif, level = lf$Motif)

ggplot(sum_stat, aes(x = ID, y = Count, fill = Motif)) +
  geom_bar(stat = "identity", position = "fill", width = 0.4) +
  scale_fill_manual(values = rev(brewer.pal(8, "Set2"))) +
  labs(x = "", y = "Percentage") +
  theme_classic() +
  mytheme +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45))
