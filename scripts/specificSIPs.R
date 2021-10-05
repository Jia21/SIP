#----------------------------------------------------------------------------------------------#
# Figure SX: gene expression of specific SIPs (Neutrophil specific used as S6)
# creates 
# requires *_specific_genes_allInfo.txt, GTEx_mat_tissue_median
#----------------------------------------------------------------------------------------------#
library(data.table)
library(tidyr)
library(ggplot2)

setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/data/")
# setwd("/proj/yunligrp/users/lagler/SIP/data/")

getData <- function(celltype){
  specific <- fread(paste0(celltype, "_specific_genes_allInfo.txt"), drop=c(2,4:6))
  specificlong <- melt(specific, id.vars = 1:2, 
                     variable.name = "tissue", value.name="TPM")
  colnames(specificlong)[[2]] <- "geneID"; colnames(specificlong)[[1]] <- "geneName"
  # only keep genes with expression info
  specificlong <- specificlong[is.na(TPM)==F]
  
  # /proj/yunligrp/users/jwen/SIP/post/GTEx_mat_tissue_median
  gtex <- fread("GTEx_mat_tissue_median", drop=12)
  colnames(gtex)[[1]] <- "geneID"; colnames(gtex)[[2]] <- "geneName"
  gtex2 <- gtex[geneName %in% specific$baitName]
  gtexlong <- melt(gtex2, id.vars = 1:2, 
                   variable.name = "tissue", value.name="TPM")
  
  # combine data
  exprdata <- rbind(gtexlong, specificlong)
  # ~3 duplicated genes with 0s
  exprdata2 <- exprdata[, .("medTPM"=median(TPM, na.rm=T)), by=c("geneName", "tissue")]
  setDT(exprdata2)[, color := ifelse(tissue=="Whole Blood", "Whole Blood",
                                     ifelse(tissue %in% c("Ery", "MacMon", "MK",
                                                          "nCD4", "Neu"), "Blood Cell", "Other Tissue"))]
  setDT(exprdata2)[, CT := celltype]
  return(exprdata2)
}

exprdata2 <- rbind(getData("Ery"), getData("MacMon"),
                   getData("MK"), getData("nCD4"),
                   getData("Neu"))

# all specific genes
ggplot(exprdata2, aes(x=tissue, y=log10(medTPM+1), color=color)) +
  geom_violin()+
  geom_point(alpha=.3) + 
  facet_wrap(~CT, ncol=1, scales="free") +
  scale_color_manual(values=c("red", "black", "blue")) +  
  xlab("Tissue/Cell Type") + ylab("log10(TPM + 1)") +
  theme_light() + 
  theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1),
        legend.title = element_blank(), legend.position = "top", 
        plot.margin=unit(c(5.5,5.5,5.5,22),"points"))
ggsave("../figures/specificSIPs.pdf", width=20, height=40, unit="in", limitsize = F)


# average by GTEx and CT
avgGeneOther <- exprdata2[color!="Blood Cell",][, .("meanExpr"=mean(medTPM)),
                                                by=c("CT", "geneName")]
setDT(avgGeneOther)[, TopPerOther := cut(meanExpr, 
                                         quantile(meanExpr, seq(.1, 1, .1)), 
                                         labels = seq(90,10,-10),
                                         include.lowest = F), by="CT"]
# average by blood cell and CT
avgGeneBlood <- exprdata2[color=="Blood Cell",][, .("meanExpr"=mean(medTPM)),
                                                by=c("CT", "geneName")]
setDT(avgGeneBlood)[, TopPerBlood := cut(meanExpr, 
                                         quantile(meanExpr, seq(.1, 1, .1), na.rm=T), 
                                         labels = seq(90,10,-10),
                                         include.lowest = F), by="CT"]
# add back to expression
exprdata3 <- merge(exprdata2, avgGeneOther[,c(1,2,4)], by="geneName")
exprdata3 <- merge(exprdata3, avgGeneBlood[,c(1,2,4)], by="geneName")
exprdata3$TopPerOther <- factor(exprdata3$TopPerOther, levels=seq(10,90,10))
exprdata3$TopPerBlood <- factor(exprdata3$TopPerBlood, levels=seq(10,90,10))

#----------------------------------------------------------------------------------------------#
# by GTEx expression
ggplot(exprdata3[TopPerOther==20], aes(x=tissue, y=log10(medTPM+1), color=color)) +
  geom_violin()+
  geom_point(alpha=.6, size=1) +
  facet_wrap(~CT, ncol=1, scales="free") + 
  scale_color_manual(values=c("red", "black", "blue")) +  
  xlab("Tissue/Cell Type") + ylab("log10(TPM + 1)") +
  theme_light() + 
  theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1),
        legend.title = element_blank(), legend.position = "top", 
        plot.margin=unit(c(5.5,5.5,5.5,22),"points"))
ggsave("../figures/specificSIPs_meanOther_20.pdf",
       width=15, height=35, unit="in", limitsize = F)

# by Blood expression
ggplot(exprdata3[TopPerBlood==10], aes(x=tissue, y=log10(medTPM+1), color=color)) +
  geom_violin()+
  geom_point(alpha=.6, size=1) +
  facet_wrap(~CT, ncol=1, scales="free") + 
  scale_color_manual(values=c("red", "black", "blue")) +  
  xlab("Tissue/Cell Type") + ylab("log10(TPM + 1)") +
  theme_light() + 
  theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1),
        legend.title = element_blank(), legend.position = "top", 
        plot.margin=unit(c(5.5,5.5,5.5,22),"points"))
ggsave("../figures/specificSIPs_meanBlood_10.pdf",
       width=15, height=35, unit="in", limitsize = F)
#----------------------------------------------------------------------------------------------#
# Figures for Supplemental: only top 10 and 20 for one cell type

mkplot <- function(plotdata){
  ggplot(plotdata, aes(x=tissue, y=log10(medTPM+1), color=color)) +
    geom_violin()+
    geom_point(alpha=.6, size=.75) +
    scale_color_manual(values=c("red", "black", "blue")) +  
    xlab("Tissue/Cell Type") + ylab("log10(TPM + 1)") +
    theme_light() + 
    theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1),
          legend.title = element_blank(), legend.position = "top", 
          plot.margin=unit(c(5.5,5.5,5.5,22),"points"))
}  

# Neutrophil specific 
ggarrange(
  mkplot(exprdata3[TopPerBlood==10 & CT=="Neu"])+ggtitle("Top 10% of Neutrophils")+
    scale_y_continuous(breaks=seq(0,5,1), limits=c(0,5.5)),
  mkplot(exprdata3[TopPerBlood==20 & CT=="Neu"])+ggtitle("Top 20% of Neutrophils")+
    scale_y_continuous(breaks=seq(0,5,1), limits=c(0,5.5)),
  nrow=2, labels=letters[1:2], align="hv", common.legend = T,
  legend = "top"
)
ggsave("../figures/specificSIPs_topNeu.png", width=10, height=11, unit="in")

# Neutrophil specific - other tissues
ggarrange(
  mkplot(exprdata3[TopPerOther==10 & CT=="Neu"])+ggtitle("Top 10% of Other Tissues")+
    scale_y_continuous(breaks=seq(0,4,1), limits=c(0,4)),
  mkplot(exprdata3[TopPerOther==20 & CT=="Neu"])+ggtitle("Top 20% of Other Tissues")+
    scale_y_continuous(breaks=seq(0,4,1), limits=c(0,4)),
  nrow=2, labels=letters[3:4], align="hv", common.legend = T,
  legend = "top"
)
ggsave("../figures/specificSIPs_topNeuOther.png", width=10, height=11, unit="in")


# MacMon specific 
ggarrange(
  mkplot(exprdata3[TopPerBlood==10 & CT=="MacMon"])+ggtitle("Top 10% of Macrophages/Monocytes")+
    scale_y_continuous(breaks=seq(0,5,1), limits=c(0,5.5)),
  mkplot(exprdata3[TopPerBlood==20 & CT=="MacMon"])+ggtitle("Top 20% of Macrophages/Monocytes")+
    scale_y_continuous(breaks=seq(0,5,1), limits=c(0,5.5)),
  nrow=2, labels=letters[1:2], align="hv", common.legend = T,
  legend = "top"
)
ggsave("../figures/specificSIPs_topMacMon.png", width=10, height=11, unit="in")

# MacMon specific - other tissues
ggarrange(
  mkplot(exprdata3[TopPerOther==10 & CT=="MacMon"])+ggtitle("Top 10% of Other Tissues")+
    scale_y_continuous(breaks=seq(0,4,1), limits=c(0,4)),
  mkplot(exprdata3[TopPerOther==20 & CT=="MacMon"])+ggtitle("Top 20% of Other Tissues")+
    scale_y_continuous(breaks=seq(0,4,1), limits=c(0,4)),
  nrow=2, labels=letters[3:4], align="hv", common.legend = T,
  legend = "top"
)
ggsave("../figures/specificSIPs_topMacMonOther.png", width=10, height=11, unit="in")

#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#