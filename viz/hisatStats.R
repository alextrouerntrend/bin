# Plant Computational Genomics Lab
# Author: Alex Trouern-Trend
# Date: May 25th, 2020

# Hisat2 output parser & stats

library(reshape2)
library(ggplot2)

setwd("~/Documents/scientist/projects/HBEF/Transcriptomics/genome-guided/hisat/")

outfa <- readLines("hisat-1099105.err")
outsp <- readLines("hisat-1099110.err")

lbfa <- c("FaFgAl1", "FaFgAl2", "FaFgAl3", "FaFgCa1", "FaFgCa2", "FaFgCo1", "FaFgCo2", "FaFgCo3" )
lbsp <- c("SpFgAl1", "SpFgAl2", "SpFgCa1", "SpFgCa2", "SpFgCa3", "SpFgCo1", "SpFgCo2", "SpFgCo3" )

hisat2_parser <- function(outfile, libnames) {
  header <- c("total_reads", "total_paired", "aligned concordantly 0 times", "aligned concordantly exactly 1 time", "aligned concordantly >1 times", "overall alignment rate")
  top <- grep("reads; of these:", outfile)
  a <- as.numeric(sub("^(\\d+).*", "\\1", outfile[top]))
  df <- data.frame("readcount" = a, stringsAsFactors = FALSE)
  nu <- 2:15
  for (i in nu) { df[i] <- as.numeric(trimws(sub("^(\\D*\\d+\\(?.\\d+).*", "\\1", outfile[top+(i-1)]))) }
  kp <- c(1,2,3,4,5,15)
  df <- df[kp]
  colnames(df) <- header
  rownames(df) <- libnames
  return(df)
}

sp <- hisat2_parser(outsp, lbsp)
fa <- hisat2_parser(outfa, lbfa)
bo <- rbind(fa, sp)

# Melt that bo fo a stack bar
bo["Library"] <- rownames(bo)
kp <-c(3,4,5,7)
bo2 <- bo[kp]
m_bo <- melt(bo2)

colors <- c("#2a4d69", "#4b86b4", "#adcbe3")

# Stacked Total
ggplot(m_bo, aes(fill=variable, y=value, x=Library)) + 
  geom_bar(position="stack", stat="identity") +
  theme_minimal() +
  scale_fill_manual(values=rev(colors)) +
  ylab("Count") +
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold", size = 10, angle = 50, hjust = 1, vjust = 1.5),
        legend.title = element_blank()) +
  ggtitle("Read Alignments Fagus grandifolia")

# Stacked Fill
ggplot(m_bo, aes(fill=variable, y=value, x=Library)) + 
  geom_bar(position="fill", stat="identity") +
  theme_minimal() +
  scale_fill_manual(values=rev(colors)) +
  ylab("Count") +
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold", size = 10, angle = 50, hjust = 1, vjust = 1.5),
        legend.title = element_blank()) +
  ggtitle("Read Alignments Fagus grandifolia")
