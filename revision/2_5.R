#----------------------------------------------------------------------------------------------#
# Revised from /proj/yunligrp/users/lagler/SIP/R/ATACseqPlot.R
# Figure 3a-b: SIPs and ATAC-seq Peak Regions
# requires SIP_*.txt, fraser.rda, and ATAC-seq peak files
# /proj/yunligrp/users/lagler/SIP/figures/ATACseqPlot.png is the table for Figure 3a.
#----------------------------------------------------------------------------------------------#
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(kableExtra)

setwd("/proj/yunligrp/users/jwen/SIP/revise/2_5")
load("/proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[31] <- "MacMon"

atac_ery <- fread("/proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/Ery_peaks.narrowPeak.gz", select=c(1:3))
atac_mon <- fread("/proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/Mono_peaks.narrowPeak.gz", select=c(1:3))
atac_ncd4 <- fread("/proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/CD4_peaks.narrowPeak.gz", select=c(1:3))
atac_mk <- fread("/proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/Mega_peaks.narrowPeak.gz", select=c(1:3))

getData <- function (cell, seqdat){
  # select significant interactions
  hic <- fraser %>% filter(!!as.name(cell) >= 5) %>% data.table
  hic <- hic[,1:11] # drop expression data
  
  # sip data
  sip <- fread(paste0("/proj/yunligrp/users/lagler/SIP/data/SIP_", cell, ".txt"))
  
  # merge with hic data to get other ends and add SIP indicator
  sipdata <- merge(hic, sip, by="baitID", allow.cartesian = T)
  # add chr to match ATAC-seq data
  setDT(sipdata)[, chr2 := paste0("chr", baitChr)]
  
  # format ATAC-seq data
  colnames(seqdat) <- c("chr", "start", "end")
  seqdat$size <- seqdat$end - seqdat$start+1
  seqdat$row <- seq(1:nrow(seqdat))
  setkey(seqdat, chr, start, end)
  
  overlap <- foverlaps(sipdata, seqdat, by.x=c("chr2", "oeStart", "oeEnd"),
                       type="any", mult="first", which=T)
  
  # add to sip_chr data
  sipdata$atacRow <- overlap # row number of ATAC-seq overlap for size merging
  sipdata$atac <- as.numeric(is.na(overlap)==F) # recode to 0/1 (false/true) overlap
  sipdata <- merge(sipdata, seqdat[, c("row", "size")],
                   by.x="atacRow", by.y="row", all.x=T, all.y=F)
  
  sipdata$chr2 <- NULL; sipdata$atacRow <- NULL
  return(as.data.table(sipdata))
  
}
#----------------------------------------------------------------------------------#
# 2x2 contigency tables of overlaps - SIP vs non-SIP - bait level
# calls overlap function
#----------------------------------------------------------------------------------#
getCountsbait <- function(cell, seqdat){
  CT2 <- getData(cell, seqdat)
  ovlps <- CT2[, .("nInteractions"=.N,
                   "nOverlaps"=sum(atac)), by=baitID]
  ovlps <- unique(merge(ovlps, CT2[, c("baitID", "SIP")], by="baitID"))
  # add cell type indicator
  setDT(ovlps)[, CT := cell]
  
  # contingency table
 # yy <- ovlps[SIP==1 & nOverlaps > 0]$nOverlaps/ovlps[SIP==1 & nOverlaps > 0]$nInteractions
  yy <- ovlps[SIP==1]$nOverlaps/ovlps[SIP==1]$nInteractions
  mean_sip_inter <- sum(ovlps[SIP==1]$nInteractions)/nrow(ovlps[SIP==1])

  #ny <- ovlps[SIP==0 & nOverlaps > 0]$nOverlaps/ovlps[SIP==0 & nOverlaps > 0]$nInteractions
  ny <- ovlps[SIP==0]$nOverlaps/ovlps[SIP==0]$nInteractions
  mean_nonsip_inter <- sum(ovlps[SIP==1]$nInteractions)/nrow(ovlps[SIP==0])
  # 
      
  # median/mean # overlaps per group
  t1 <- summary(yy)
  t2 <- summary(ny)
  counts <- cbind(c(paste0(cell,"_SIP PIR ratio"),paste0(cell,"_non_SIP PIR ratio")),rbind(c(t1,mean_sip_inter),c(t2,mean_nonsip_inter)))
  counts <- data.table(counts)
  # wilcoxon test for # overlaps
  w <- wilcox.test(yy, ny,alternative="less")$p.value
  setDT(counts)[, wilcox.pval.onesided := w]
  
  # t test for # overlaps
  t <- t.test(yy, ny,alternative="less")$p.value
  setDT(counts)[, t.pval.onesided := t]
  
  return(list("data"=list(yy,ny),
              "table"=counts))
}

ratios <- data.table(rbind(getCountsbait("Ery", atac_ery)$table,
                           getCountsbait("MacMon", atac_mon)$table,
                           getCountsbait("MK", atac_mk)$table,
                           getCountsbait("nCD4", atac_ncd4)$table))
colnames(ratios)[8] <- "Mean of the SIP/non-SIP PIR number"

ovlpData <- list(getCountsbait("Ery", atac_ery)$data,
                  getCountsbait("MacMon", atac_mon)$data,
                  getCountsbait("MK", atac_mk)$data,
                  getCountsbait("nCD4", atac_ncd4)$data)

dat1 <- rbind(rbind(cbind(ovlpData[[1]][[1]],"Ery","SIP"),cbind(ovlpData[[1]][[2]],"Ery","non-SIP")),
              rbind(cbind(ovlpData[[2]][[1]],"MacMon","SIP"),cbind(ovlpData[[2]][[2]],"MacMon","non-SIP")),
              rbind(cbind(ovlpData[[3]][[1]],"MK","SIP"),cbind(ovlpData[[3]][[2]],"MK","non-SIP")),
              rbind(cbind(ovlpData[[4]][[1]],"nCD4","SIP"),cbind(ovlpData[[4]][[2]],"nCD4","non-SIP")))

dat1 <- data.frame(dat1)
dat1$X1 <- as.numeric(as.character(dat1$X1))
mean <- data.frame(cbind(ratios$Mean,c("Ery","Ery","MacMon","MacMon","MK","MK","nCD4","nCD4"),c("SIP PIRs","no-SIP PIRs", "SIP PIRs","no-SIP PIRs","SIP PIRs","no-SIP PIRs","SIP PIRs","no-SIP PIRs")))
median <- data.frame(cbind(ratios$Median,c("Ery","Ery","MacMon","MacMon","MK","MK","nCD4","nCD4"),c("SIP PIRs","no-SIP PIRs", "SIP PIRs","no-SIP PIRs","SIP PIRs","no-SIP PIRs","SIP PIRs","no-SIP PIRs")))
mean$X1 <- as.numeric(as.character(mean$X1))
median$X1 <- as.numeric(as.character(median$X1))
dat1 <- data.table(dat1)
    plotd <- ggplot(dat1, aes(x=X3, y = X1, fill=X3)) +
          geom_violin() + scale_fill_brewer(palette="Paired") + ylab("% of PIRs overlapping ATAC-seq peaks") + 
          xlab("") +
          geom_point(data=dat1[, .("median"=median(X1)), by=X3], aes(y=median), show.legend = F,size = 3)+ 
          geom_point(data=dat1[, .("mean"=mean(X1)), by=X3], aes(y=mean), show.legend = F,size = 3,shape = 17)+ facet_wrap(.~X2) +
          theme_bw() + 
          theme(axis.text=element_text(size=20),axis.title=element_text(size=20)) +
          theme(legend.title = element_blank(),legend.text = element_text(size=20))+
          theme(strip.text.x = element_text(size = 15))
ggsave("ATAC.png", height = 9, width=12, units="in")
ggsave("ATAC.pdf", height =9, width=12, units="in")

write.table(ratios,"pir_atac_seq_ratio_comp_includingNonOlap",quote =F,sep="\t",col.names = T,row.names = F)