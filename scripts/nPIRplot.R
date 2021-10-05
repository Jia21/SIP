#----------------------------------------------------------------------------------------------#
# Figure 2(a-b): violin plots of number of PIRs and CHiCAGO scores
# requires fraser.rda and SIP_*.txt (from defineSIPlot.R)
#----------------------------------------------------------------------------------------------#
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/data/")
#setwd("/proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/") #fraser.rda
#----------------------------------------------------------------------------------------------#
# Find # of PIRs per bait and CHiCAGO scores
#----------------------------------------------------------------------------------------------#
# formatted for intra-chromosomal interactions and distance < 2MB
load("fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[31] <- "MacMon"

getData <- function (cell){
  # select significant interactions
  hic <- fraser %>% filter(!!as.name(cell) >= 5) %>% data.table
  # add chicago score column 
  hic$chicago <- hic[ , eval(as.name(cell))]
  
  # count number of PIRs by bait
  # and mean/median CHiCAGO score per bait
  bybait <- hic[, .("N"=.N,
                    "avgC"=mean(chicago),
                    "medC"=median(chicago)), by="baitID"]
  
  # merge hic data with SIP identifier
  # /proj/yunligrp/users/lagler/SIP/data/
  sip <- fread(paste0("SIP_", cell, ".txt"))
  bybait2 <- merge(bybait, sip, by="baitID")
  
  # add cell type identifier
  setDT(bybait2)[, CT := cell]
  
  return(bybait2)
}

# number of significant interactions per bait
# mean/median chicago score per bait
bybait <- rbind(getData("Ery"), getData("MacMon"), getData("MK"),
                getData("nCD4"), getData("Neu"))

# significant differences in each case (p < 2.2e-16)
wilcox.test(bybait[CT=="Neu" & SIP==1, N], bybait[CT=="Neu" & SIP==0, N])
wilcox.test(bybait[CT=="Neu" & SIP==1, medC], bybait[CT=="Neu" & SIP==0, medC])
# median estimates
bybait[, .("medianN"=median(N),
           "medianMedC"=median(medC)), by=c("CT", "SIP")]
#----------------------------------------------------------------------------------------------#
# Make violin plots of distributions
#----------------------------------------------------------------------------------------------#

# number of significant interactions per bait
ninter <- ggplot(data=bybait, aes(x=CT, y=log10(N), fill=as.factor(SIP))) + 
  geom_violin()+
  geom_point(data=bybait[, .("median"=median(N)), by=c("CT", "SIP")],
             aes(y=log10(median)), position =position_dodge(.9), show.legend = F)+
  scale_fill_brewer(palette="Paired", direction=1, labels = c("non-SIP", "SIP"))+ 
  ylab("Number of Interactions per Bait (log10)") + 
  xlab("Cell Type Group") + 
  theme_bw()+theme(legend.title = element_blank())

# median chicago score per bait
medChicago <- ggplot(data=bybait, aes(x=CT, y=log10(medC), fill=as.factor(SIP))) + 
  geom_violin(trim=T) +
  geom_point(data=bybait[, .("median"=median(medC)), by=c("CT", "SIP")],
             aes(y=log10(median)), position =position_dodge(.9), show.legend = F)+
  scale_fill_brewer(palette="Paired", direction=1, labels = c("non-SIP", "SIP"))+ 
  ylab("Median CHiCAGO Score per Bait (log10)") + 
  xlab("Cell Type Group") + 
  theme_bw()+theme(legend.title = element_blank())

ggarrange(ninter, medChicago, align="hv",
          labels=letters[1:2], common.legend = T)
ggsave("../figures/violinSIP.pdf", width=9, height=4)
ggsave("../figures/violinSIP.png", width=9, height=4)

#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#