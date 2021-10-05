#----------------------------------------------------------------------------------------------#
# Figure 3a-b: SIPs and ATAC-seq Peak Regions
# requires SIP_*.txt, fraser.rda, and ATAC-seq peak files
#----------------------------------------------------------------------------------------------#
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(kableExtra)
setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/data/")
#----------------------------------------------------------------------------------------------#
# read in data 
#----------------------------------------------------------------------------------------------#

# pcHiC data
# formatted for intra-chromosomal interactions and distance < 2MB
load("fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[31] <- "MacMon"

# ATAC-seq data
setwd("/Users/tmlagler/OneDrive/Lab/SIP/ATACseq/")
# setwd("/proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/")
atac_ery <- fread("Ery_peaks.narrowPeak.gz", select=c(1:3))
atac_mon <- fread("Mono_peaks.narrowPeak.gz", select=c(1:3))
atac_ncd4 <- fread("CD4_peaks.narrowPeak.gz", select=c(1:3))
atac_mk <- fread("Mega_peaks.narrowPeak.gz", select=c(1:3))

#----------------------------------------------------------------------------------------------#
# find overlaps with ATAC-seq by CT
#----------------------------------------------------------------------------------------------#
setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/data/")
# setwd("/proj/yunligrp/users/lagler/SIP/data/")
getData <- function (cell, seqdat){
  # select significant interactions
  hic <- fraser %>% filter(!!as.name(cell) >= 5) %>% data.table
  hic <- hic[,1:11] # drop expression data
  
  # sip data
  sip <- fread(paste0("SIP_", cell, ".txt"))
  
  # merge with hic data to get other ends and add SIP indicator
  sipdata <- merge(hic, sip, by="baitID", allow.cartesian = T)
  # add chr to match ATAC-seq data
  setDT(sipdata)[, chr2 := paste0("chr", baitChr)]
  
  # format ATAC-seq data
  colnames(seqdat) <- c("chr", "start", "end")
  seqdat$size <- seqdat$end - seqdat$start
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
  yy <- nrow(ovlps[SIP==1 & nOverlaps > 0])
  yn <- nrow(ovlps[SIP==1 & nOverlaps == 0])
  ny <- nrow(ovlps[SIP==0 & nOverlaps > 0])
  nn <- nrow(ovlps[SIP==0 & nOverlaps == 0])
  
  # 2x2 table of counts plus ratio
  counts <- data.table("CT" = rep(cell, 2),
                       "type" = c("SIP", "non-SIP"),
                       "ATACseqOvlp" = c(yy, ny),
                       "noATACseqOvlp" = c(yn, nn))
  setDT(counts)[, Ratio := ATACseqOvlp/(ATACseqOvlp+noATACseqOvlp)]
  
  # chisq test
  p <- chisq.test(counts[,3:4])$p.value
  setDT(counts)[, chsq.pval := p]
  
  # median/mean # overlaps per group
  medOvlp <- ovlps[, .("med" = median(nOverlaps),
                       "avg" = mean(nOverlaps)),, by="SIP"]
  counts$medOvlp <- medOvlp[c(2,1), med]
  counts$avgOvlp <- medOvlp[c(2,1), avg]
  
  # wilcoxon test for # overlaps
  w <- wilcox.test(ovlps[SIP==1, nOverlaps], ovlps[SIP==0, nOverlaps])$p.value
  setDT(counts)[, wilcox.pval := w]
  
  # t test for # overlaps
  t <- t.test(ovlps[SIP==1, nOverlaps], ovlps[SIP==0, nOverlaps])$p.value
  setDT(counts)[, t.pval := t]
  
  return(list("data"=ovlps,
              "table"=counts))
}


#----------------------------------------------------------------------------------------------#
# call function for each cell type. make plots
#----------------------------------------------------------------------------------------------#

# ratio of ATAC-seq overlaps (SIP v non-SIP)
# ratio of any overlap
ratios <- data.table(rbind(getCountsbait("Ery", atac_ery)$table,
                           getCountsbait("MacMon", atac_mon)$table,
                           getCountsbait("MK", atac_mk)$table,
                           getCountsbait("nCD4", atac_ncd4)$table))

# used to add significance marks
maxRatio <- ratios[, .("max"=max(Ratio)), by="CT"]
maxRatio$type <- "SIP"

# flip plot to combine with GWAS PIR
ratios$CT <- factor(ratios$CT, levels=c("Neu", "nCD4", "MK", "MacMon", "Ery"))
propPlot <- 
  ggplot(data=ratios, aes(y=Ratio, x=CT, fill=type))+
    coord_flip(ylim=c(.5,1.03))+
  geom_bar(stat="identity", position=position_dodge())+
  geom_point(data=maxRatio, aes(y=max+0.03), shape=8, color="red", show.legend=F)+
  scale_fill_brewer(palette="Paired", direction=1, labels = c("non-SIP", "SIP"))+ 
  geom_label(aes(label=round(Ratio,3), group=factor(type)),
            vjust=0.3, hjust=1.25,
            position=position_dodge(width=0.9),
            color="black", fill="white", size=2.5) +
  ylab("Proportion Overlapping ATAC-seq") + 
  xlab("Cell Type Group") +
  theme_bw() + theme(legend.title = element_blank())

# number of ATAC-seq overlaps
ovlpData <- rbind(getCountsbait("Ery", atac_ery)$data,
                  getCountsbait("MacMon", atac_mon)$data,
                  getCountsbait("MK", atac_mk)$data,
                  getCountsbait("nCD4", atac_ncd4)$data)

# used to add significance marks
maxN <- ovlpData[, .("max"=max(nOverlaps)), by="CT"]
maxN$SIP <- 1

# flip plot to combine with GWAS PIR
ovlpData$CT <- factor(ovlpData$CT, levels=c("Neu", "nCD4", "MK", "MacMon", "Ery"))

nOvlpPlot <- 
  ggplot(data=ovlpData, aes(y=nOverlaps, x=CT, fill=as.factor(SIP)))+
    coord_flip(ylim=c(0,max(ovlpData$nOverlaps)))+
  geom_boxplot(outlier.alpha=.6)+
  geom_point(data=maxN, aes(y=max+1), shape=8, color="red", show.legend=F)+
  scale_fill_brewer(palette="Paired", direction=1, labels = c("non-SIP", "SIP"))+ 
  ylab("Number of ATAC-seq Overlaps") + 
  xlab("Cell Type Group") + 
  theme_bw()+theme(legend.title = element_blank())

# combine plots
ggarrange(propPlot, nOvlpPlot, nrow=1, common.legend = T, labels = c("a", "b"))
ggsave("./figures/ATAC.png", height = 4, width=9, units="in")
ggsave("./figures/ATAC.pdf", height = 4, width=9, units="in")

#----------------------------------------------------------------------------------------------#
# make supplementary table
#----------------------------------------------------------------------------------------------#

# keep formatting separate incase want to change digits
ratiosK <- ratios
ratiosK$chsq.pval <- format(ratiosK$chsq.pval, digits=2) # format p-values
ratiosK$wilcox.pval <- format(ratiosK$wilcox.pval, digits=2)
ratiosK$t.pval <- format(ratiosK$t.pval, digits=2)
ratiosK$chsq.pval[c(2,4,6,8)] <- NA # remove every other value
ratiosK$wilcox.pval[c(2,4,6,8)] <- NA
ratiosK$t.pval[c(2,4,6,8)] <- NA
ratiosK$Ratio <- format(ratiosK$Ratio, digits=2, drop0trailing = F)
ratiosK$avgOvlp <- format(ratiosK$avgOvlp, digits=2, drop0trailing = F)
ratiosK[,1] <- NA

options(knitr.kable.NA = "")

# wilcoxon p-values all 0
ratiosK$wilcox.pval <- NULL
kable(ratiosK, col.names = c("", "Bait Type", "# Overlap",
                             "# No Overlap", "Ratio", "Chi-sq p-value",
                          "Med. # Ovlp.", "Avg. # Ovlp.", "t p-value"),
      booktabs=T, "latex", align="llcccccccc", 
      format.args = list(big.mark=",")) %>%
  pack_rows(index = c("Ery" = 2, "MacMon" = 2, "MK" = 2, "nCD4"=2)) %>%
  kable_styling(latex_options = c("scale_down")) %>%
  save_kable("./figures/ATACseqTable.pdf")

#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#
