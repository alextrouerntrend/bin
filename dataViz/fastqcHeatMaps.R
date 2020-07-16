# Author: Alex Trouern-Trend
# Date: May 18th, 2020

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Documents/scientist/projects/HBEF/Transcriptomics/read_qc/")

sumtr <- read.csv("summ_trimmed.tsv", header = FALSE, sep = "\t")
sumog <- read.csv("summ_original.tsv", header = FALSE, sep = "\t")

sumtr_stat <- dcast(sumtr, V2 ~ V1)
sumog_stat <- dcast(sumog, V2 ~ V1)


tidy_sumtr <- sumtr %>%
spread(V2, V1)

tidy_sumog <- sumog %>%
  spread(V2, V1)

springas <- sumtr[grep("SpAs", sumtr$V3),]
springfg <- sumtr[grep("SpFg", sumtr$V3),]
fallas <- sumtr[grep("FaAs", sumtr$V3),]
fallfg <- sumtr[grep("FaFg", sumtr$V3),]

colors <- c("#d1495b", "#66a182", "#edae49")

# Heatmap for Spring Acer saccharum
ggplot(springas, aes(x = V2, y = V3, fill = factor(V1))) + 
  geom_tile() + 
  scale_fill_manual(values=colors) +
  theme(text=element_text(family="mono"),
        axis.text.x = element_text(angle = 50, hjust = 1),
        axis.title = element_blank(),
        legend.title = element_blank(),
        panel.background = element_blank(),
        title = element_text(face = "bold")) +
  ggtitle("Spring: Acer saccharum")

# Heatmap for Spring Fagus grandifolia
ggplot(springfg, aes(x = V2, y = V3, fill = factor(V1))) + 
  geom_tile() + 
  scale_fill_manual(values=colors) +
  theme(text=element_text(family="mono"),
        axis.text.x = element_text(angle = 50, hjust = 1),
        axis.title = element_blank(),
        legend.title = element_blank(),
        panel.background = element_blank(),
        title = element_text(face = "bold")) +
  ggtitle("Spring: Fagus grandifolia")

# Heatmap for Fall Acer saccharum
ggplot(fallas, aes(x = V2, y = V3, fill = factor(V1))) + 
  geom_tile() + 
  scale_fill_manual(values=colors) +
  theme(text=element_text(family="mono"),
        axis.text.x = element_text(angle = 50, hjust = 1),
        axis.title = element_blank(),
        legend.title = element_blank(),
        panel.background = element_blank(),
        title = element_text(face = "bold")) +
  ggtitle("Fall: Acer saccharum")

# Heatmap for Fall Fagus grandifolia
ggplot(fallfg, aes(x = V2, y = V3, fill = factor(V1))) + 
  geom_tile() + 
  scale_fill_manual(values=colors) +
  theme(text=element_text(family="mono"),
        axis.text.x = element_text(angle = 50, hjust = 1),
        axis.title = element_blank(),
        legend.title = element_blank(),
        panel.background = element_blank(),
        title = element_text(face = "bold")) +
  ggtitle("Fall: Fagus grandifolia")

