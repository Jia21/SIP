#----------------------------------------------------------------------------------------------#
# Figure 5: Dot Plot of Fold Enrichment
# requires SIP_*.txt and allSIPs_geneInfo.txt
# requires *_specific_genes_allInfo.txt for CT specific SIPs
#----------------------------------------------------------------------------------------------#
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(coin) # independence_test()
library(scales)
setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/data/")
#setwd("/proj/yunligrp/users/lagler/SIP/data/") 
#----------------------------------------------------------------------------------------------#
# Get fold enrichment and trend test results for each cell type
#----------------------------------------------------------------------------------------------#

# expression data with baitNames from (getSIPgeneInfo.R)
expr <- fread("allSIPs_geneInfo.txt")
expr <- expr[ENSEMBL_GENEID!="" & is.na(Ery)==F]

getData <- function(cell, specific=F){
  # add entropy data to SIP call data
  if(specific==T){sipfile <- paste0("SIP_", cell, "_specific.txt")}else{
    sipfile <- paste0("SIP_", cell, ".txt")
    specfile <- fread(paste0(cell, "_specific_genes_allInfo.txt"))
  }
  sip <- fread(sipfile)
  sip$rank <- NULL
  
  sip2 <- merge(sip, expr[,c("baitID", "baitName", "Ery", "MacMon","MK",
                             "nCD4", "Neu")], by="baitID")
  
  if(specific==T){
    setDT(sip2)[, SIP := ifelse(baitName %in% specfile$baitName, 1, 0)]
    sip2[,4] <- NULL
    # sip$SIP <- sip[,4] # use specific indicator as SIP indicator
    # sip[,4] <- NULL
    }
  # wide to long
  siplong <- melt(sip2, value.vars = c("MacMon", "MK", "nCD4", "Ery", "Neu"), 
                  id.vars = c("baitID", "intscore", "SIP", "baitName"),
                  value.name = "expression", variable.name = "CT")
  # rank of expression values by gene
  setDT(siplong)[, rank := frank(-expression), by=c("baitID", "baitName")]
  
  # trend test
  trend <- independence_test(SIP ~ rank, data=siplong[CT==cell, ], alternative="less")
  
  # get fold enrichment and p-value
  foldData <- data.table("rank"=1:5)
  for(r in 1:5){
    rankR <- siplong[rank==r & CT== cell]
    foldData$nGenes[r] <- nrow(rankR) # total number of genes at rank r
    foldData$nSIP[r] <- nrow(rankR[SIP==1]) # total number of SIPs at rank r
    foldData$prop[r] <- foldData$nSIP[r]/foldData$nGenes[r] # proportion of SIPs at rank r
  }
  for(r in 1:5){
    # fold enrichment relative to 5th rank
    foldData$FE[r] <- foldData$prop[r]/foldData$prop[5]
    # chi-sq p-value
    foldData$chisq.pval[r] <- prop.test(c(foldData[r, nSIP], foldData[5, nSIP]),
                                        c(foldData[r, nGenes], foldData[5, nGenes]))$p.value
  }
  
  # add trend test p-value to output
  setDT(foldData)[, trend.pval := coin::pvalue(trend)]
  
  # add cell type
  setDT(foldData)[, CT := cell]
  
  return(foldData)
}

#----------------------------------------------------------------------------------------------#
# call function for each cell type and make plot - all SIPs
#----------------------------------------------------------------------------------------------#

enrich <- rbind(getData("Ery"), getData("MacMon"), getData("MK"),
                getData("nCD4"), getData("Neu"))
#enrich$log10chip <- -log10(enrich$chisq.pval)
# reorder cell types so Ery is on top
enrich$CT <- factor(enrich$CT, levels=c("Neu", "nCD4", "MK", "MacMon", "Ery"))

# format trend test p-values
trend <- unique(enrich[, c("CT", "trend.pval")])
setDT(trend)[, fp := format(trend.pval, digits=2, drop0trailing = F)]

# make plot
ggplot(data=enrich, aes(x=rank, y=CT, size=FE, fill=chisq.pval)) + 
  geom_point(color="black", pch=21) + 
  scale_fill_gradient2(low="yellow", mid="yellow", high="skyblue3",
                       breaks=pretty_breaks(n=4)(min(enrich$chisq.pval):max(enrich$chisq.pval)),
                       limits=c(min(enrich$chisq.pval), max(enrich$chisq.pval)),
                       midpoint=0.05, labs(title="Chi-sq p-value")) +
  scale_size_continuous(breaks=c(1, 1.4, 1.8),
                        labs(title="Fold Enrichment"))+
  xlab("Gene Expression Ranking") + ylab("")+
  scale_x_continuous(labels=c("1st", "2nd", "3rd", "4th", "5th", "Trend Test\np-value"),
                     breaks=c(1, 2, 3, 4, 5, 6), position="top")+
  geom_vline(xintercept = 5.5, color="black", linetype="dashed")+
  theme_minimal() +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 1, 0, 0, "cm"),
        #panel.border = element_rect(color="black", fill=NA),
        legend.position = "bottom")+
  # Trend test p-values
  annotate("text", x=6, y=trend$CT,
           label=trend$fp)

ggsave("../figures/foldEnrichPlot.png", width=7, height=4)
ggsave("../figures/foldEnrichPlot.pdf", width=7, height=4)

#----------------------------------------------------------------------------------------------#
# call function for each cell type and make plot - CT specific SIPs
# all SIPs preferred 
#----------------------------------------------------------------------------------------------#

enrich_S <- rbind(getData("Ery", specific=T), getData("MacMon", specific=T),
                getData("MK", specific=T), getData("nCD4", specific=T),
                getData("Neu", specific=T))
#enrich$log10chip <- -log10(enrich$chisq.pval)
# reorder cell types so Ery is on top
enrich_S$CT <- factor(enrich_S$CT, levels=c("Neu", "nCD4", "MK", "MacMon", "Ery"))

# format trend test p-values
trend_S <- unique(enrich_S[, c("CT", "trend.pval")])
setDT(trend_S)[, fp := format(trend.pval, digits=2, drop0trailing = F)]

# make plot
ggplot(data=enrich_S, aes(x=rank, y=CT, size=FE, fill=chisq.pval)) + 
  geom_point(color="black", pch=21) + 
  scale_fill_gradient2(low="yellow", mid="yellow", high="skyblue3",
                       #breaks=pretty_breaks(n=4)(min(enrich_S$chisq.pval):max(enrich_S$chisq.pval)),
                       limits=c(min(enrich_S$chisq.pval), max(enrich_S$chisq.pval)),
                       midpoint=0.05, labs(title="Chi-sq p-value")) +
  scale_size_continuous(breaks=c(1, 1.5, 2),
                        labs(title="Fold Enrichment"))+
  xlab("Gene Expression Ranking") + ylab("")+
  scale_x_continuous(labels=c("1st", "2nd", "3rd", "4th", "5th", "Trend Test\np-value"),
                     breaks=c(1, 2, 3, 4, 5, 6), position="top")+
  geom_vline(xintercept = 5.5, color="black", linetype="dashed")+
  theme_minimal() +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 1, 0, 0, "cm"),
        #panel.border = element_rect(color="black", fill=NA),
        legend.position = "bottom")+
  # Trend test p-values
  annotate("text", x=6, y=trend_S$CT,
           label=trend_S$fp)

ggsave("../figures/foldEnrichPlotSpecific.png", width=7, height=4)
ggsave("../figures/foldEnrichPlotSpecific.pdf", width=7, height=4)

#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#
