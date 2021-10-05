################################################################################
### Using the output from /proj/yunligrp/users/lagler/SIP/R/getSIPgeneInfo.R ###
###                   By Jia Wen 09/2021                               #########
################################################################################
rm(list=ls())
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
setwd("/proj/yunligrp/users/jwen/SIP/revise")
allSIPgeneInfo <- read.table("/proj/yunligrp/users/lagler/SIP/data/allSIPs_geneInfo.txt",fill=T,sep="\t",header =T)
# tpm <- function(vec){
#      vec <- na.omit(vec)
#      tmp_vec <- vec/sum(vec) * 1e6
#      return(tmp_vec)
# }
# MK <- tpm(allSIPgeneInfo$MK)
# Ery <- tpm(allSIPgeneInfo$Ery)
# MacMon <- tpm(allSIPgeneInfo$Mon)
# nCD4 <- tpm(allSIPgeneInfo$nCD4)
# Neu <- tpm(allSIPgeneInfo$Neu)
cell <- c("Ery","MacMon","MK","nCD4","Neu")

t <- NULL
for (c in 4:8){
    cat(colnames(allSIPgeneInfo)[c],"\t")
    sipgene <- allSIPgeneInfo[,cell[c-3]][allSIPgeneInfo[,c]==1]
    nonsipgene <- allSIPgeneInfo[,cell[c-3]][allSIPgeneInfo[,c]==0]
    t1 <- data.frame(sipgene,paste0(cell[c-3],".SIP"), cell[c-3], "SIP")
    colnames(t1) <- c("value","type","CT","category")
    t2 <- data.frame(nonsipgene,paste0(cell[c-3],".non.SIP"), cell[c-3], "non-SIP")
    colnames(t2) <- c("value","type","CT","category")
    t <- rbind(t,t1,t2)
}
t <- data.frame(t)
t_noNA <- na.omit(t)

#----------------------------------------------------------------------------------------------#
# Make violin plots of distributions
#----------------------------------------------------------------------------------------------#
t_noNA <- data.table(t_noNA)
t_noNA$category <- factor(t_noNA$category,levels= c("non-SIP", "SIP"))
ggplot(data=t_noNA, aes(x=CT, y=log(value), fill=category)) + 
  geom_violin(trim=T) +
  geom_point(data=t_noNA[, .("median"=median(log(value))), by=c("CT", "category")], aes(y=median), position =position_dodge(.9), show.legend = F)+
  scale_fill_brewer(palette="Paired", direction=1) + 
  ylab("log(Gene expression)") + 
  xlab("Cell Type Group") + 
  theme_bw()+theme(legend.title = element_blank())
ggsave("./1_7/Figure6b_violinSIP_nonSIP_exp.png", width=9, height=7)
ggsave("./1_7/Figure6b_violinSIP_nonSIP_exp.pdf", width=9, height=9)

test <- NULL
for (c in 1:length(cell)){
    p <- wilcox.test(as.numeric(as.character(t_noNA[t_noNA$CT==cell[c] & t_noNA$category=="SIP", "value"])), 
                     as.numeric(as.character(t_noNA[t_noNA$CT==cell[c] & t_noNA$category=="non-SIP", "value"])),alternativ="greater")$p.value
    test <- rbind(test,c(cell[c],p))
}
#      [,1]     [,2]                  
# [1,] "Ery"    "1.16001146806389e-56"
# [2,] "MacMon" "7.13820422438885e-13"
# [3,] "MK"     "9.10771106902095e-28"
# [4,] "nCD4"   "2.95294684977504e-88"
# [5,] "Neu"    "5.59926596638126e-86"