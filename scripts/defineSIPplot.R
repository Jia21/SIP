#----------------------------------------------------------------------------------------------#
# Figure 1: Hockey stick plots of Interaction Scores (defining SIPs)
# requires fraser.rda and annotgenes.txt
#----------------------------------------------------------------------------------------------#
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)
setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/data/")
#setwd("/proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/") #fraser.rda
#----------------------------------------------------------------------------------------------#
# pcHiC data and selected genes for annotation
#----------------------------------------------------------------------------------------------#
# formatted for intra-chromosomal interactions and distance < 2MB
load("fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[31] <- "MacMon"

# selected cell type group-specific genes
# /proj/yunligrp/users/lagler/SIP/data/annotgenes.txt
genes <- fread("annotgenes.txt", drop="rank")

#----------------------------------------------------------------------------------------------#
# Calculate Cumulative Interaction Score for each Cell Type Group
# find cutoff for defining SIPs
# plot curve. annotate selected CT-specific genes
#----------------------------------------------------------------------------------------------#

mkPlot <- function(cell, cellname){
  # filter for significant interactions based on CT
  hic <- fraser %>% filter(!!as.name(cell) >= 5) %>% data.table
  
  # sum CHiCAGO scores by bait
  hic$ct <- hic[ , eval(as.name(cell))]
  score <- hic[, .(intscore = sum(ct)), by="baitID"]
  # rank cumulative interaction scores
  setorder(score, intscore)
  score$rank <- seq(1:nrow(score))
  
  # find inflection point (i.e. define SIPs)
  # normalize rank and score
  setDT(score)[, norm_rank := rank/max(rank)]
  setDT(score)[, norm_score := intscore/max(intscore)]
  setDT(score)[, rot1 := 1/sqrt(2)*norm_rank + 1/sqrt(2)*norm_score]
  setDT(score)[, rot2 := -1/sqrt(2)*norm_rank + 1/sqrt(2)*norm_score]
  
  RefPoint <- score[rot2==min(rot2), rank]
  RefValue <- score[rot2==min(rot2), intscore] # any score > is a SIP
  # add SIP indicator
  setDT(score)[, SIP := ifelse(intscore > RefValue, 1, 0)]
  NumSIP <- sum(score$SIP)
  
  # # add rank
  annot <- merge(genes[CT==cell],
                 score[, c("baitID", "rank")], by="baitID")
  setorder(annot, -intscore)
  
  
  # hockey stick plot
  hplot <- ggplot(score, aes(x=rank, y=intscore)) + 
    geom_point(size=1) + 
    geom_vline(xintercept=RefPoint, color="blue", linetype="dashed") +
    
    # annotate genes
    geom_point(data=annot, color="red", size=1.5) +
    geom_text_repel(data=annot, aes(label=baitName), segment.color = "darkgray", 
                    point.padding=1, box.padding=2, nudge_y = 400, nudge_x=-3000) + 
    
    xlab("Rank of Promoter Bait") + ylab("Cumulative Interaction Score") + 
    scale_x_continuous(breaks=c(0, 5000, 10000))+
    ggtitle(cellname)+
    theme_bw()
  
  # return plot, number of SIPs found, and data
  return(list("plot" = hplot, "nSIPs" = c(cellname, NumSIP),
              "data" = score[, c("baitID", "intscore", "rank", "SIP")]))
}

#----------------------------------------------------------------------------------------------#
# Call Function for each Cell Type Group
#----------------------------------------------------------------------------------------------#
ery <- mkPlot("Ery", "Erythrocyte")
macmon <- mkPlot("MacMon", "Macrophage/Monocyte")
mk <- mkPlot("MK", "Megakaryocyte")
ncd4 <- mkPlot("nCD4", "Naive CD4 T-cell")
neu <- mkPlot("Neu", "Neutrophil")

# make table of SIPs per CT
nSIPs <- rbind(ery$nSIPs, macmon$nSIPs, mk$nSIPs, ncd4$nSIPs, neu$nSIPs)
nSIPs[,2] <- format(as.numeric(nSIPs[,2]), big.mark = ",")
colnames(nSIPs) <- c("Cell Type", "# of SIPs")
ggtable <- ggtexttable(nSIPs, theme=ttheme("lBlack", base_size=12))

# arrange all eCDF plots and table
ggarrange(ery$plot, macmon$plot, mk$plot, ncd4$plot, neu$plot, 
          ggtable,
          nrow=2, ncol=3, labels=letters[1:6],
          common.legend = T)
ggsave("../figures/Figure1_hockeySIP.pdf", width=11, height=7)
ggsave("../figures/Figure1_hockeySIP.png", width=11, height=7)

# save SIP data
# /proj/yunligrp/users/lagler/SIP/data/
fwrite(ery$data, "SIP_Ery.txt", col.names=T, sep="\t")
fwrite(macmon$data, "SIP_MacMon.txt", col.names=T, sep="\t")
fwrite(mk$data, "SIP_MK.txt", col.names=T, sep="\t")
fwrite(ncd4$data, "SIP_nCD4.txt", col.names=T, sep="\t")
fwrite(neu$data, "SIP_Neu.txt", col.names=T, sep="\t")

#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#
